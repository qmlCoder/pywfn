from pywfn.base import Mol,Atom
import numpy as np
from pywfn.data import Elements
from functools import lru_cache
elements=Elements()
from pywfn.utils import printer
from pywfn.atomprop import lutils,AtomCaler
from pywfn import maths
from typing import Literal
from pywfn.maths import CM2PM

Chrgs=Literal['mulliken','lowdin','hirshfeld']

class Calculator(AtomCaler):
    def __init__(self,mol:"Mol"):
        self.logTip:str=''
        self.mol=mol
        self.chrg:Chrgs='mulliken'
    
    def calculate(self,chrg:Chrgs,PM=None)->np.ndarray:
        if chrg=='mulliken':
            return self.mulliken(PM)
        if chrg=='lowdin':
            return self.lowdin(PM)
        if chrg=='hirshfeld':
            return self.hirshfeld()
    
    def mulliken(self,num:bool=False,PM:np.ndarray=None):
        """
        计算目录mulliken电荷
        num：是否只保留电子数
        """
        # 计算密度矩阵
        self.logTip='Mulliken电荷分布'
        if PM is None:PM=self.mol.PM
        # 矩阵乘法的迹的加和=矩阵对应元素乘积之和
        PS=PM@self.mol.SM
        EV=np.diagonal(PS)
        atoms=self.mol.atoms
        charges=np.zeros(len(atoms))
        for a,atom in enumerate(atoms):
            a1,a2=atom.obtBorder
            elect=EV[a1:a2].sum()
            if num:
                charges[a]=elect
            else:
                charges[a]=atom.atomic-elect
        return charges
    
    def lowdin(self,num:bool=False,PM:np.ndarray=None):
        """
        计算每个原子的lowdin电荷
        """
        self.logTip='Lowdin电荷分布'
        if PM is None:PM=self.mol.PM
        SM=self.mol.SM
        # 计算矩阵的1/2次方
        v, Q=np.linalg.eig(SM)
        V=np.diag(v)
        V_=V**0.5
        Q_=np.linalg.inv(Q)
        SM_half=Q@(V_@Q_)
        SPS=SM_half@(PM@SM_half)
        eleNums=np.diagonal(SPS)
        charges=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            eleNum=eleNums[u:l].sum()
            charges[a]=eleNum if num else atom.atomic-eleNum
        return charges
    
    def sapce(self,num:bool=False,PM:np.ndarray=None):
        """计算空间电荷,DFT格点数值积分"""
        from pywfn.spaceProp import density
        if PM is None:PM=self.mol.PM
        caler=density.Calculator(self.mol)
        charges=np.zeros(len(self.mol.atoms))
        for a,atom in enumerate(self.mol.atoms):
            dens=caler.atmDens(atom.idx)
            eleNum=np.sum(dens)
            charges[a]=eleNum if num else atom.atomic-eleNum
        return charges
    
    def hirshfeld(self):
        """
        计算原子的Hirshfeld电荷
        """
        from pywfn.spaceProp import density
        from pywfn.data import radDens
        denCaler=density.Calculator(self.mol)
        molPos,molWei=denCaler.molPos
        
        npos=len(molPos)
        natm=len(self.mol.atoms)
        pdens=np.zeros(npos) # 前体电子密度
        fdensl=[]
        for atom in self.mol.atoms:
            radius=np.linalg.norm(molPos-atom.coord,axis=1)
            fdens=radDens.get_radDens(atom.atomic,radius)
            pdens+=fdens
            fdensl.append(fdens)
        Ps=np.zeros(natm)
        Zs=np.zeros(natm)
        mdens=denCaler.molDens_atm(molPos)
        for a,atom in enumerate(self.mol.atoms): # 计算每一个原子的电荷
            radius=np.linalg.norm(molPos-atom.coord,axis=1)
            fdens=fdensl[a]
            Wa=np.divide(fdens,pdens,out=np.zeros_like(fdens),where=pdens!=0)
            Zs[a]=np.sum(fdens*molWei) # 自由态下核电荷数
            Ps[a]=np.sum(Wa*mdens*molWei) # 真实体系的电子布局
        chargs=Zs-Ps
        
        return chargs
    
    
    def dirCharge(self,chrg:Chrgs,atms:list[int],dirs:list[np.ndarray]=None)->np.ndarray:
        """计算不同方向的电荷[n,5](atm,x,y,z,val)"""
        fatms,fdirs=fit_dirs(self.mol,atms,dirs)
        assert len(fatms)==len(fdirs),"长度需要一致"
        dirVal=np.zeros(shape=(len(fdirs),5))
        obts=self.mol.O_obts
        for d in range(len(fdirs)):
            fatm=fatms[d]
            fdir=fdirs[d]
            CMp=self.mol.projCM(obts,[fatm],[fdir],False,False) # 获取投影后的轨道系数，单个原子投影到指定方向
            PMp=CM2PM(CMp,obts,self.mol.oE)
            if chrg=='mulliken':
                val=self.mulliken(num=True,PM=PMp)[fatm-1]
            elif chrg=='lowdin':
                val=self.lowdin(num=True,PM=PMp)[fatm-1]
            x,y,z=fdir
            dirVal[d]=[fatm,x,y,z,val]
        return dirVal
    
    def piElectron(self,chrg:Chrgs='mulliken'):
        """计算π电子"""
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        atms=[]
        dirs=[]
        idxs=[]
        for idx,atom in enumerate(self.mol.atoms):
            normal=dirCaler.normal(atom.idx) # 原子的法向量
            if normal is None:continue
            atms.append(atom.idx)
            dirs.append(normal)
            idxs.append(idx)
        CMp=self.mol.projCM(self.mol.O_obts,atms,dirs,False,False) # 所有能投影的原子同时投影各自的法向量
        PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
        if chrg=='mulliken':
            val=self.mulliken(num=True,PM=PMp)
        elif chrg=='lowdin':
            val=self.lowdin(num=True,PM=PMp)
        atms=np.array(atms).reshape(-1,1)
        dirs=np.vstack(dirs)
        val=val[idxs].reshape(-1,1)
        return np.hstack([atms,dirs,val])
    
    def onShell(self):
        from pywfn.utils import parse_intList
        chrgMap={'':'mulliken','1':'mulliken','2':'lowdin'}
        while True:
            printer.options('原子电荷',{
                '1':'mulliken电荷',
                '2':'lowdin电荷',
                '3':'方向性电荷',
                '4':'π电子分布',
            })
            opt=input('请选择电荷类型:')
            if opt=='1':
                charges=self.mulliken()
                for i,v in enumerate(charges):
                    print(f'{i+1}: {v}')
            elif opt=='2':
                charges=self.lowdin()
                for i,v in enumerate(charges):
                    print(f'{i+1}: {v}')
            elif opt=='3':
                print('1. Mulliken[*]; 2. lowdin')
                opt=input('选择电荷类型:')
                if opt not in chrgMap.keys():return
                chrg=chrgMap[opt]
                numStr=input('输入要计算的原子：')
                atms=parse_intList(numStr)
                charges=self.dirCharge(chrg,atms)
                for a,x,y,z,v in charges:
                    print(f'{int(a):>3d}({x:>6.2f},{y:>6.2f},{z:>6.2f}):{v:>8.4f}')
            elif opt=='4':
                print('1. Mulliken[*]; 2. lowdin')
                opt=input('选择电荷类型:')
                if opt not in chrgMap.keys():return
                chrg=chrgMap[opt]
                result=self.piElectron(chrg)
                for i,val in enumerate(result):
                    print(f'{i+1:>3d}:{val:>8.4f}')
            else:
                break
    
def fit_dirs(mol:Mol,atms:list[int],dirs:list[np.ndarray]):
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