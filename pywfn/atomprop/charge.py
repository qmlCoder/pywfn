from pywfn.base import Mole,Atom
import numpy as np
from pywfn.data.elements import elements
from functools import lru_cache
from pywfn.utils import printer
from pywfn.gridprop import density,dftgrid
from pywfn import maths
from typing import Literal
from pywfn.maths import CM2PM
from pywfn.maths.mol import projCM
from pywfn.cli import Shell
from collections import defaultdict
from pywfn.atomprop import lutils

Chrgs=Literal['mulliken','lowdin','space','hirshfeld']

class Calculator():
    def __init__(self,mol:"Mole"):
        self.logTip:str=''
        self.mol=mol
        self.form:str='charge' # 输出格式，电子数或电荷数 number|charge
        self.PM=self.mol.PM.copy() # 计算时使用的密度矩阵
    
    def spin(self,ctype:str='mulliken')->np.ndarray:
        """计算所有原子的自旋"""
        if not self.mol.open: #闭壳层自旋肯定为0
            return np.zeros(self.mol.atoms.natm)
        nmat=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        occs=[int(e) for e in self.mol.obtOccs] # 记录原本的占据情况
        CMa=self.mol.CM[:,:nmat//2].copy() # alpha轨道系数
        CMb=self.mol.CM[:,nmat//2:].copy() # beta轨道系数
        
        PMa=CM2PM(CMa,occs[:nmat//2],1) # alpha密度矩阵
        PMb=CM2PM(CMb,occs[nmat//2:],1) # beta密度矩阵
        SM=self.mol.SM # 重叠矩阵
        atmuls=self.mol.atoms.atmuls # 原子的轨道数量
        match ctype:
            case 'mulliken':
                electA=lutils.mulliken(PMa,SM,atmuls)
                electB=lutils.mulliken(PMb,SM,atmuls)
                return electA-electB
            case _:
                raise ValueError(f'未知的点和类型{ctype}')
    
    def mulliken(self,form:str='charge')->np.ndarray:
        """
        计算mulliken电荷
        """
        elects=lutils.mulliken(self.mol.PM,self.mol.SM,self.mol.atoms.atmuls)
        match form:
            case 'number':
                return elects
            case 'charge':
                return np.array(self.mol.atoms.atomics)-elects
            case _:
                raise ValueError(f'未知的输出格式{form}')
    
    def lowdin(self,form:str='charge')->np.ndarray:
        """
        计算每个原子的lowdin电荷
        """
        elects=lutils.lowdin(self.mol.PM,self.mol.SM,self.mol.atoms.atmuls)
        atmics=np.array(self.mol.atoms.atomics)
        match form:
            case 'number':
                return elects
            case 'charge':
                return atmics-elects
            case _:
                raise ValueError(f'未知的输出格式{form}')
    
    def space(self,form:str='charge')->np.ndarray:
        """计算空间电荷,DFT格点数值积分"""
        gridCaler=dftgrid.Calculator(self.mol)
        grids,weits=gridCaler.molGrid()
        densCaler=density.Calculator(self.mol)
        densCaler.PM=self.PM # 设为当前PM
        elects=np.zeros(len(self.mol.atoms))
        atmics=np.array(self.mol.atoms.atomics)
        for a,atom in enumerate(self.mol.atoms):
            dens=densCaler.atmDens(grids,[atom.idx])
            eleNum=np.sum(dens*weits)
            elects[a]=eleNum
        match form:
            case 'number':
                return elects
            case 'charge':
                return atmics-elects
            case _:
                raise ValueError(f'未知的输出格式{form}')
    
    def hirshfeld(self,form:str='charge')->np.ndarray:
        """计算原子的Hirshfeld电荷"""
        from pywfn.data import radDens
        from pywfn.gridprop import dftgrid
        from pywfn.gridprop import density
        gridCaler=dftgrid.Calculator(self.mol) # 格点计算器
        gridCaler.nrad=75
        gridCaler.nsph=434
        densCaler=density.Calculator(self.mol) # 电子密度计算器
        natm=self.mol.atoms.natm
        chargs=np.zeros(natm)
        atmics=np.array(self.mol.atoms.atomics)
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
            mdens=densCaler.molDens(atmGrids,0)[0]
            res=np.divide(idens,pdens,where=pdens!=0)*mdens*weits*gcuts
            chargs[i]=np.sum(res)
        match form:
            case 'number':
                return chargs
            case 'charge':
                return atmics-chargs
            case _:
                raise ValueError(f'未知的输出格式{form}')

    def dirElects(self,atms:list[int],dirs:np.ndarray,ctype:str)->np.ndarray:
        """计算不同方向的电子数"""
        assert len(atms)==len(dirs),"原子与方向数量要相同"
        obts=self.mol.O_obts
        self.numForm=True
        CMp=projCM(self.mol,obts,atms,dirs,False,False) # 获取投影后的轨道系数，单个原子投影到指定方向
        PMp=CM2PM(CMp,obts,self.mol.oE)
        match ctype:
            case 'mulliken':
                vals=lutils.mulliken(PMp,self.mol.SM,self.mol.atoms.atmuls)
            case 'lowdin':
                vals=lutils.lowdin(PMp,self.mol.SM,self.mol.atoms.atmuls)
            case _:
                raise ValueError(f'未知电荷类型{ctype}')
        return vals
    
    def piElects(self,ctype:str='mulliken'):
        """
        计算π电子数
        每个原子的方向为`法向量`的`方向电子数`即为`π电子数`
        """
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol)
        atms=[]
        dirs=[]
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx)
            if normal is None:continue
            atms.append(atom.idx)
            x,y,z=normal
            dirs.append([x,y,z])
        dirs=np.array(dirs)
        return atms,dirs,self.dirElects(atms,dirs,ctype)
    
    def piDecomElects(self,ctype='mulliken'): # 使用轨道分解方法计算pi电子分布，可以包含D轨道
        from pywfn.orbtprop import decom

        CMt=decom.Calculator(self.mol).pi_decom('atom')
        PMt=CM2PM(CMt,self.mol.O_obts,self.mol.oE) # 变换的密度矩阵
        match ctype:
            case 'mulliken':
                return lutils.mulliken(PMt,self.mol.SM,self.mol.atoms.atmuls)
            case 'lowdin':
                return lutils.lowdin(PMt,self.mol.SM,self.mol.atoms.atmuls)
            case _:
                raise ValueError(f'未知电荷类型{ctype}')

    
def fit_dirs(mol:Mole,atms:list[int],dirs:list[np.ndarray]|None):
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