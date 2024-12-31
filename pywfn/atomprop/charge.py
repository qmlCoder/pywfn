from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data.elements import elements
from functools import lru_cache
from pywfn.utils import printer
from pywfn.spaceprop import density,dftgrid
from pywfn import maths
from typing import Literal
from pywfn.maths import CM2PM
from pywfn.maths.mol import projCM
from pywfn.shell import Shell
from collections import defaultdict

Chrgs=Literal['mulliken','lowdin','space','hirshfeld']

class Calculator():
    def __init__(self,mol:"Mol"):
        self.logTip:str=''
        self.mol=mol
        self.numForm:bool=False # 输出格式，电子数或电荷数 number|charge
        self.PM=self.mol.PM.copy() # 计算时使用的密度矩阵
        
    def charge(self,ctype:str)->np.ndarray:
        """计算四种基础电荷之一

        Args:
            chrg (Chrgs): 电荷类型

        Returns:
            np.ndarray: 原子电荷
        """
        if ctype=='mulliken':
            return self.mulliken()
        elif ctype=='lowdin':
            return self.lowdin()
        elif ctype=='space':
            return self.sapce()
        elif ctype=='hirshfeld':
            return self.hirshfeld()
        else:
            raise ValueError('unknown charge type')
    
    def mulliken(self)->np.ndarray:
        """
        计算目录mulliken电荷
        num：是否只保留电子数
        """
        # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
        PS=self.PM@self.mol.SM
        EV=np.diagonal(PS) # 矩阵的对角元素
        atoms=self.mol.atoms
        elects=np.zeros(len(atoms))
        for a,atom in enumerate(atoms):
            u,l=atom.obtBorder
            elect=EV[u:l].sum()
            elects[a]=elect
            # charges[a]=atom.atomic-elect
        if self.numForm:
            return elects
        else:
            return np.array(self.mol.atoms.atomics)-elects
    
    def lowdin(self)->np.ndarray:
        """
        计算每个原子的lowdin电荷
        """
        PM=self.PM
        SM=self.mol.SM
        # 计算矩阵的1/2次方
        v, Q=np.linalg.eig(SM)
        V=np.diag(v)
        V_=V**0.5
        Q_=np.linalg.inv(Q)
        SM_half=Q@(V_@Q_)
        SPS=SM_half@(PM@SM_half)
        eleNums=np.diagonal(SPS)
        elects=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            eleNum=eleNums[u:l].sum()
            elects[a]=eleNum
        if self.numForm:
            return elects
        else:
            return np.array(self.mol.atoms.atomics)-elects
    
    def sapce(self)->np.ndarray:
        """计算空间电荷,DFT格点数值积分"""
        gridCaler=dftgrid.Calculator(self.mol)
        grid,weit=gridCaler.molGrid()
        densCaler=density.Calculator(self.mol)
        # PMo=densCaler.PM # 备份旧的PM
        densCaler.PM=self.PM # 设为当前PM
        elects=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            dens=densCaler.atmDens(grid,[atom.idx])
            eleNum=np.sum(dens*weit)
            elects[a]=eleNum
        # densCaler.PM=PMo # 恢复旧的PM
        if self.numForm:
            return elects
        else:
            return np.array(self.mol.atoms.atomics)-elects
    
    def hirshfeld(self)->np.ndarray:
        """计算原子的Hirshfeld电荷"""
        from pywfn.data import radDens
        from pywfn.spaceprop import dftgrid
        from pywfn.spaceprop import density
        gridCaler=dftgrid.Calculator(self.mol) # 格点计算器
        densCaler=density.Calculator(self.mol) # 电子密度计算器
        natm=self.mol.atoms.natm
        Ps=np.zeros(natm)
        Zs=np.zeros(natm)
        chargs=np.zeros(natm)
        for i,atom in enumerate(self.mol.atoms):
            grid,weit=gridCaler.atmGrid(atom.idx)
            npos=len(grid)
            pdens=np.zeros(npos) # 前体电子密度
            fdensl=[] # 原子自由电子密度
            for j,atom in enumerate(self.mol.atoms):
                radius=np.linalg.norm(grid-atom.coord,axis=1) # 所有格点到当前原子的距离
                fdens=radDens.get_radDens(atom.atomic,radius) # 自由原子电子密度
                pdens+=fdens
                fdensl.append(fdens)
            mdens,_,_=densCaler.molDens(grid,0) #分子电子密度
            res=fdensl[i]/pdens*mdens*weit
            chargs[i]=np.sum(res)
        if self.numForm:
            return chargs
        else:
            return np.array(self.mol.atoms.atomics)-chargs

    def dirElectron(self,atms:list[int],dirs:np.ndarray,ctype:str)->np.ndarray:
        """计算不同方向的电子[n,5](atm,x,y,z,val)"""
        assert len(atms)==len(dirs),"原子与方向数量要相同"
        obts=self.mol.O_obts
        PMo=self.PM.copy() # 记录旧的密度矩阵
        self.numForm=True
        CMp=projCM(self.mol,obts,atms,dirs,False,False) # 获取投影后的轨道系数，单个原子投影到指定方向
        PMp=CM2PM(CMp,obts,self.mol.oE)
        self.PM=PMp # 修改密度矩阵
        vals=self.charge(ctype)
        self.PM=PMo # 恢复旧的密度矩阵
        return vals
    
    def piElectron(self,ctype:str):
        """计算π电子"""
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        atmDirs={} # 记录每一个原子的方向
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx) # 原子的法向量
            if normal is None:continue # 没有法向量就跳过
            atmDirs[atom.idx]=normal
        atms=list(atmDirs.keys())
        dirs=np.array([e for e in atmDirs.values()])
        self.numForm=True
        CMp=projCM(self.mol,self.mol.O_obts,atms,dirs,False,False) # 所有能投影的原子同时投影各自的法向量
        PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
        PMo=self.PM.copy() # 备份老的密度矩阵
        self.PM=PMp # 将密度矩阵替换为投影后的密度矩阵
        eleNums=self.charge(ctype)
        # print(eleNums)
        result=[]
        for atom in self.mol.atoms:
            atm=atom.idx
            if atm in atms:
                x,y,z=atmDirs[atm]
            else:
                x,y,z=[0,0,0]
            ele=eleNums[atm-1]
            result.append([atm,x,y,z,ele])
        self.PM=PMo
        return np.array(result)
    
    def piElectDecom(self): # 使用轨道分解方法计算pi电子分布，可以包含D轨道
        from pywfn.orbtprop import decom

        CMt=decom.Calculator(self.mol).pi_decom('atom')
        PMt=CM2PM(CMt,self.mol.O_obts,self.mol.oE) # 变换的密度矩阵
        self.PM=PMt
        self.numForm=True
        return self.mulliken()

    
    def onShell(self,shell:Shell):
        from pywfn.utils import parse_intList
        chrgMap={'':'mulliken','1':'mulliken','2':'lowdin','3':'space','4':'hirshfeld'}
        chrgStr='1. Mulliken[*]; 2. lowdin; 3. sapce; 4. hirshfeld'
        while True:
            printer.options('原子电荷',{
                '1':'mulliken电荷',
                '2':'lowdin电荷',

                '3':'空间积分电荷',
                '4':'hirshfeld电荷',

                '5':'方向电子分布',
                '6':'π 电子分布',
            })
            opt=input('请选择计算类型: ')
            if opt=='1': # Mulliken
                charges=self.mulliken()
                for i,val in enumerate(charges):
                    print(f'{i+1:>3d}: {val:>8.4f}')
                print(f'sum:{np.sum(charges)}')

            elif opt=='2': # Lowdin
                charges=self.lowdin()
                for i,val in enumerate(charges):
                    print(f'{i+1:>3d}: {val:>8.4f}')
                print(f'sum:{np.sum(charges)}')
            
            elif opt=='3': # 空间积分电荷
                charges=self.sapce()
                for i,val in enumerate(charges):
                    print(f'{i+1:>3d}: {val:>8.4f}')
                print(f'sum:{np.sum(charges)}')

            elif opt=='4': # hirshfeld电荷
                charges=self.hirshfeld()
                for i,val in enumerate(charges):
                    print(f'{i+1:>3d}: {val:>8.4f}')
                print(f'sum:{np.sum(charges)}')

            elif opt=='5': # 方向电子
                print(chrgStr)
                opt=input('选择电荷类型: ')
                if opt not in chrgMap.keys():return
                chrg=chrgMap[opt]
                # numStr=input('输入原子编号: ')
                atm=shell.input.Integ(tip='输入原子编号: ',count=1)[0]
                dir=shell.input.Float(tip='输入原子向量: ',count=3)
                assert dir is not None,'原子向量输入不正确'
                dirs=np.array(dir).reshape(1,3)
                elects=self.dirElectron(atms=[atm],dirs=dirs,ctype=chrg)
                ele=elects[atm-1]
                x,y,z=dir
                print(f'{atm:>3d} ({x:>6.2f},{y:>6.2f},{z:>6.2f}):{ele:>8.4f}')
            
            elif opt=='6': # pi电子
                print(chrgStr)
                opt=input('选择电荷类型: ')
                if opt not in chrgMap.keys():return
                chrg=chrgMap[opt]
                elects=self.piElectron(chrg)
                for i,(idx,x,y,z,val) in enumerate(elects):
                    print(f'{i+1:>3d}:{val:>8.4f}')
                print(f'sum:{elects[:,-1].sum()}')

            else:
                break
    
def fit_dirs(mol:Mol,atms:list[int],dirs:list[np.ndarray]|None):
    """
    矫正方向，如果没有指定方向的话，计算每个原子可能的反应方向
    """
    if dirs is None:
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(mol)
        fatms=[]
        fdirs=[]
        for atm in atms:
            resDirs=dirCaler.reactions(atm)
            fdirs.append(resDirs)
            fatms+=[atm]*len(resDirs)
        fdirs=np.vstack(fdirs)
    else:
        fatms=atms
        fdirs=dirs
    return fatms,fdirs