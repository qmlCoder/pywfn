from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data.elements import elements
from functools import lru_cache
from pywfn.utils import printer
from pywfn.atomprop import lutils,AtomCaler
from pywfn.spaceProp import density,dftgrid
from pywfn import maths
from typing import Literal
from pywfn.maths import CM2PM
from pywfn.maths.mol import projCM

Chrgs=Literal['mulliken','lowdin','space','hirshfeld']

class Calculator(AtomCaler):
    def __init__(self,mol:"Mol"):
        self.logTip:str=''
        self.mol=mol
        self.chrg:Chrgs='mulliken'
        self.form:str='charge' # 输出格式，电子数或电荷数
        self.PM=self.mol.PM.copy() # 计算时使用的密度矩阵
        
    def charge(self,chrg:Chrgs)->np.ndarray:
        """计算四种基础电荷之一

        Args:
            chrg (Chrgs): 电荷类型

        Returns:
            np.ndarray: 原子电荷
        """
        if chrg=='mulliken':
            return self.mulliken()
        elif chrg=='lowdin':
            return self.lowdin()
        elif chrg=='space':
            return self.sapce()
        elif chrg=='hirshfeld':
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
        if self.form=='number':
            return elects
        if self.form=='charge':
            return np.array(self.mol.atoms.atomics)-elects
        # return charges
    
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
        if self.form=='number':
            return elects
        if self.form=='charge':
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
            dens=densCaler.atmDens(atom.idx,grid)
            eleNum=np.sum(dens*weit)
            elects[a]=eleNum
        # densCaler.PM=PMo # 恢复旧的PM
        if self.form=='number':
            return elects
        if self.form=='charge':
            return np.array(self.mol.atoms.atomics)-elects
    
    def hirshfeld(self)->np.ndarray:
        """
        计算原子的Hirshfeld电荷
        """
        from pywfn.spaceProp import density
        from pywfn.data import radDens
        gridCaler=dftgrid.Calculator(self.mol)
        
        densCaler=density.Calculator(self.mol)
        grid,weit=gridCaler.molGrid()
        densCaler.set_grid(grid)
        PMo=densCaler.PM # 备份旧的PM
        densCaler.PM=self.PM # 设为当前PM
        
        npos=len(grid)
        natm=self.mol.atoms.natm
        pdens=np.zeros(npos) # 前体电子密度
        fdensl=[] # 原子自由电子密度
        # print(time.time())
        for atom in self.mol.atoms:
            radius=np.linalg.norm(grid-atom.coord,axis=1)
            fdens=radDens.get_radDens(atom.atomic,radius)
            pdens+=fdens
            fdensl.append(fdens)

        Ps=np.zeros(natm)
        Zs=np.zeros(natm)
        
        mdens=densCaler.molDens_lib()

        for a,atom in enumerate(self.mol.atoms): # 计算每一个原子的电荷
            radius=np.linalg.norm(grid-atom.coord,axis=1)
            fdens=fdensl[a]
            Wa=np.divide(fdens,pdens,out=np.zeros_like(fdens),where=pdens!=0)
            Zs[a]=np.sum(fdens*weit) # 自由态下核电荷数
            Ps[a]=np.sum(Wa*mdens*weit) # 真实体系的电子布局
        chargs=Zs-Ps
        densCaler.PM=PMo # 恢复旧的PM
        if self.form=='number':
            return np.array(self.mol.atoms.atomics)-chargs
        if self.form=='charge':
            return chargs
    
    def dirElectron(self,atms:list[int],dirs:list[np.ndarray],ctype:str)->np.ndarray:
        """计算不同方向的电子[n,5](atm,x,y,z,val)"""
        # fatms,fdirs=fit_dirs(self.mol,atms,dirs) # 矫正原子索引和方向，使数量相等
        assert len(atms)==len(dirs),"原子与方向数量要相同"
        # result=np.zeros(shape=(len(dirs),5))
        obts=self.mol.O_obts
        PMo=self.PM.copy() # 记录旧的密度矩阵
        self.form='number'

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
        dirs=list(atmDirs.values())
        
        CMp=projCM(self.mol,self.mol.O_obts,atms,dirs,False,False) # 所有能投影的原子同时投影各自的法向量
        PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
        PMo=self.PM.copy() # 备份老的密度矩阵
        self.PM=PMp # 将密度矩阵替换为投影后的密度矩阵
        self.form='number'
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
    
    def onShell(self):
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
                numStr=input('输入原子编号: ')
                atms=parse_intList(numStr,start=1)
                elects=self.dirElectron(chrg,atms)
                for a,x,y,z,v in elects:
                    print(f'{int(a):>3d}({x:>6.2f},{y:>6.2f},{z:>6.2f}):{v:>8.4f}')
                # print(f'sum:{elects[:,-1].sum()}')
            
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
            resDirs=dirCaler.reaction(atm)
            fdirs.append(resDirs)
            fatms+=[atm]*len(resDirs)
        fdirs=np.vstack(fdirs)
    else:
        fatms=atms
        fdirs=dirs
    return fatms,fdirs