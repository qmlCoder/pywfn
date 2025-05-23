"""
轨道分解
"""
from collections import defaultdict
import numpy as np


from pywfn.base import Mole
from pywfn.atomprop import direction

class Calculator:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol

    def pi_decom(self,dtype:str)->np.ndarray: # 分解出pi分子轨道
        dirCaler=direction.Calculator(self.mol)
        CMt=self.mol.CM.copy()
        nmat=self.mol.CM.shape[0]
        # Ts=[]

        keeps=[[0],[0,0,1],[0,0,0,1,1,1]]
        bases=dirCaler.bases()
        nobt=self.mol.CM.shape[1]
        for o in range(nobt): # 遍历轨道
            coefDict=defaultdict(list) # 系数字典
            for i in range(nmat):
                iatm=self.mol.atoAtms[i]
                ishl=self.mol.atoShls[i]
                iang=self.mol.atoAngs[i]
                key=(iatm,ishl,iang)
                coefDict[key].append(self.mol.CM[i,o])

            for key,val in coefDict.items():
                iatm,ishl,iang=key
                rcoefs=np.array(val) # 原始系数
                if iatm in bases.keys():
                    T=bases[iatm]
                    tcoefs=decomOrbitals(T,rcoefs,keeps[iang],dtype)
                else:
                    tcoefs=np.zeros_like(rcoefs)
                coefDict[key]=tcoefs.tolist() # type: ignore
            values=list(coefDict.values())
            CMt[:,o]=np.concatenate(values)
        return CMt


def decomOrbitals(T:np.ndarray,coefs:np.ndarray,keeps:list[int],dtype:str):
    match len(coefs):
        case 1:
            tcoefs = decomOrbitalS(T,coefs,keeps)
        case 3:
            tcoefs = decomOrbitalP(T,coefs,keeps)
        case 6:
            if dtype=='atom':
                tcoefs=coefs*np.array([0.,0.,0.,1.,1.,1.])
            elif dtype=='bond':
                tcoefs = decomOrbitalD(T,coefs,keeps)
            else:
                raise Exception('未知轨道类型')
        case 10:
            pass
        case _:
            # return coefs
            raise Exception('未知轨道类型')
    # print(f'{coefs}->{tcoefs}')
    return tcoefs

def decomOrbitalS(T:np.ndarray,coefs:np.ndarray,keeps:list[int]):
    if keeps[0]:
        return coefs
    else:
        return np.array([0.])

# 分解P轨道
def decomOrbitalP(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    """分解P轨道

    Args:
        T (np.ndarray): 基坐标，每一行代表一个方向
        rcoefs (np.ndarray): 原始函数空间的基函数系数
        keeps (list[int]): 保留的角动量

    Returns:
        np.ndarray: 分解之后的轨道系数
    """
    Mr=np.linalg.inv(T)
    tcoefs=Mr@rcoefs # 根据函数空间基组1下的系数获取函数空间基组2下的系数
    tcoefs*=np.array(keeps) # 根据角动量保留的系数
    Mi=np.linalg.inv(Mr)
    fcoefs=Mi@tcoefs # 根据修改后的函数空间基组2下的系数得到函数空间基组1下的系数
    return fcoefs

# 分解D轨道
def decomOrbitalD(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    M=np.array([
        [T[0,0]**2, T[0,0]*T[0,1], T[0,0]*T[0,2], T[0,1]**2, T[0,1]*T[0,2], T[0,2]**2],
        [T[1,0]**2, T[1,0]*T[1,1], T[1,0]*T[1,2], T[1,1]**2, T[1,1]*T[1,2], T[1,2]**2],
        [T[2,0]**2, T[2,0]*T[2,1], T[2,0]*T[2,2], T[2,1]**2, T[2,1]*T[2,2], T[2,2]**2],
        [2*T[0,0]*T[1,0], (T[0,0]*T[1,1]+T[0,1]*T[1,0]), (T[0,0]*T[1,2]+T[0,2]*T[1,2]), 2*T[0,1]*T[1,1], (T[0,1]*T[1,2]+T[0,2]*T[1,1]), 2*T[0,2]*T[1,2]],
        [2*T[0,0]*T[2,0], (T[0,0]*T[2,1]+T[0,1]*T[2,0]), (T[0,0]*T[2,2]+T[0,2]*T[2,0]), 2*T[0,1]*T[2,1], (T[0,1]*T[2,2]+T[0,2]*T[2,1]), 2*T[0,2]*T[2,2]],
        [2*T[1,0]*T[2,0], (T[1,0]*T[2,1]+T[1,1]*T[2,0]), (T[1,0]*T[2,2]+T[1,2]*T[2,0]), 2*T[1,1]*T[2,1], (T[1,1]*T[2,2]+T[1,2]*T[2,1]), 2*T[1,2]*T[2,2]],
    ])
    Mr=np.linalg.inv(M)
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    tcoefs*=np.array(keeps)
    fcoefs=Mi@tcoefs
    return fcoefs