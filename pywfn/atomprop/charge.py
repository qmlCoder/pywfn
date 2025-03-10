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
    
    def spin(self,chrg:str='mulliken')->np.ndarray:
        """计算所有原子的自旋"""
        if not self.mol.open: #闭壳层自旋肯定为0
            return np.zeros(self.mol.atoms.natm)
        nmat=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        occs_old=self.mol.obtOccs # 记录原本的占据情况
        
        a_occs=occs_old.copy() # 当你需要修改一个变量的时候，
        b_occs=occs_old.copy()
        a_occs[nmat:]=[False]*nmat
        b_occs[:nmat]=[False]*nmat
        elects=[]
        for occs in (a_occs,b_occs):
            self.mol._props.set('obtOccs',occs)
            elect=self.charge(chrg) # 计算电荷分布
            elects.append(elect)
        a_elect,b_elect=elects
        print(f'atm:{"Na":>10},{"Nb":>10}')
        for i,(ela,elb) in enumerate(zip(a_elect,b_elect)):
            print(f'{i+1:>3d}:{ela:>10.4f},{elb:>10.4f}')
        print('-'*25)
        # 恢复分子属性
        self.mol._props.set('obtOccs',occs_old)
        return -(a_elect-b_elect)
    
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
        gridCaler.nrad=75
        gridCaler.nsph=434
        densCaler=density.Calculator(self.mol) # 电子密度计算器
        natm=self.mol.atoms.natm
        chargs=np.zeros(natm)
        grids,weits,gcuts=gridCaler.dftGrid(1)
        npos=len(grids)
        for i,atom in enumerate(self.mol.atoms): # 循环每个原子
            atmGrids=grids+atom.coord.reshape(1,3)
            pdens=np.zeros(npos) # 前体电子密度
            idens=np.zeros(shape=(natm,npos)) # 原子自由电子密度
            for j,atom in enumerate(self.mol.atoms):
                dists=np.linalg.norm(atmGrids-atom.coord.reshape(1,3),axis=1) # 所有格点到当前原子的距离
                fdens=radDens.get_radDens(atom.atomic,dists) # 自由原子电子密度
                pdens+=fdens #*gcuts
                if i==j:idens=fdens
            mdens=densCaler.atmDens(atmGrids).sum(axis=0)
            res=np.divide(idens,pdens,where=pdens!=0)*mdens*weits*gcuts
            # for j in range(npos):
            #     x,y,z=atmGrids[j]
            chargs[i]=np.sum(res)
            # print(i,chargs[i])
            print(f'\r{i+1}/{natm}',end='')
        print('')
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