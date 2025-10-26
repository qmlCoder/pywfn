"""
轨道分解
"""
from collections import defaultdict
import numpy as np
from pprint import pprint

from pywfn.base.mole import Mole
from pywfn.atomprop import direction

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole

    def pi_decom(self,dtype:str)->np.ndarray: # 分解出pi分子轨道
        dirCaler=direction.Calculator(self.mole)
        CMt=self.mole.CM.copy()
        nmat=self.mole.CM.shape[0]
        # Ts=[]

        keeps=[[0],[0,0,1],[0,0,0,1,1,1]]
        bases=dirCaler.bases()
        # print('bases',bases)
        nobt=self.mole.CM.shape[1]
        for o in range(nobt): # 遍历轨道
            coefDict=defaultdict(list) # 系数字典
            for i in range(nmat):
                iatm=self.mole.atoAtms[i]
                ishl=self.mole.atoShls[i]
                iang=self.mole.atoAngs[i]
                key=(iatm,ishl,iang)
                coefDict[key].append(self.mole.CM[i,o].item())
            # pprint(coefDict)
            for key,val in coefDict.items():
                iatm,ishl,iang=key
                rcoefs=np.array(val) # 原始系数
                if iatm in bases.keys():
                    T=bases[iatm]
                    atom=self.mole.atom(iatm)
                    if atom.is_linear:
                        tcoefs=decomOrbitals(T,rcoefs,[0,1,1],dtype)
                    else:
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

def decomOrbitalF(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    T11=T[0,0]
    T12=T[0,1]
    T13=T[0,2]
    T21=T[1,0]
    T22=T[1,1]
    T23=T[1,2]
    T31=T[2,0]
    T32=T[2,1]
    T33=T[2,2]
    MF=np.zeros((10,10))
    MF[0][0]=T11**3
    MF[0][1]=T12**3
    MF[0][2]=T13**3
    MF[0][3]=T11*T12**2
    MF[0][4]=T11**2*T12
    MF[0][5]=T11**2*T13
    MF[0][6]=T11*T13**2
    MF[0][7]=T12*T13**2
    MF[0][8]=T12**2*T13
    MF[0][9]=T11*T12*T13
    MF[1][0]=T21**3
    MF[1][1]=T22**3
    MF[1][2]=T23**3
    MF[1][3]=T21*T22**2
    MF[1][4]=T21**2*T22
    MF[1][5]=T21**2*T23
    MF[1][6]=T21*T23**2
    MF[1][7]=T22*T23**2
    MF[1][8]=T22**2*T23
    MF[1][9]=T21*T22*T23
    MF[2][0]=T31**3
    MF[2][1]=T32**3
    MF[2][2]=T33**3
    MF[2][3]=T31*T32**2
    MF[2][4]=T31**2*T32
    MF[2][5]=T31**2*T33
    MF[2][6]=T31*T33**2
    MF[2][7]=T32*T33**2
    MF[2][8]=T32**2*T33
    MF[2][9]=T31*T32*T33
    MF[3][0]=3*T11*T21**2
    MF[3][1]=3*T12*T22**2
    MF[3][2]=3*T13*T23**2
    MF[3][3]=T11*T22**2 + 2*T12*T21*T22
    MF[3][4]=2*T11*T21*T22 + T12*T21**2
    MF[3][5]=2*T11*T21*T23 + T13*T21**2
    MF[3][6]=T11*T23**2 + 2*T13*T21*T23
    MF[3][7]=T12*T23**2 + 2*T13*T22*T23
    MF[3][8]=2*T12*T22*T23 + T13*T22**2
    MF[3][9]=T11*T22*T23 + T12*T21*T23 + T13*T21*T22
    MF[4][0]=3*T11**2*T21
    MF[4][1]=3*T12**2*T22
    MF[4][2]=3*T13**2*T23
    MF[4][3]=2*T11*T12*T22 + T12**2*T21
    MF[4][4]=T11**2*T22 + 2*T11*T12*T21
    MF[4][5]=T11**2*T23 + 2*T11*T13*T21
    MF[4][6]=2*T11*T13*T23 + T13**2*T21
    MF[4][7]=2*T12*T13*T23 + T13**2*T22
    MF[4][8]=T12**2*T23 + 2*T12*T13*T22
    MF[4][9]=T11*T12*T23 + T11*T13*T22 + T12*T13*T21
    MF[5][0]=3*T11**2*T31
    MF[5][1]=3*T12**2*T32
    MF[5][2]=3*T13**2*T33
    MF[5][3]=2*T11*T12*T32 + T12**2*T31
    MF[5][4]=T11**2*T32 + 2*T11*T12*T31
    MF[5][5]=T11**2*T33 + 2*T11*T13*T31
    MF[5][6]=2*T11*T13*T33 + T13**2*T31
    MF[5][7]=2*T12*T13*T33 + T13**2*T32
    MF[5][8]=T12**2*T33 + 2*T12*T13*T32
    MF[5][9]=T11*T12*T33 + T11*T13*T32 + T12*T13*T31
    MF[6][0]=3*T11*T31**2
    MF[6][1]=3*T12*T32**2
    MF[6][2]=3*T13*T33**2
    MF[6][3]=T11*T32**2 + 2*T12*T31*T32
    MF[6][4]=2*T11*T31*T32 + T12*T31**2
    MF[6][5]=2*T11*T31*T33 + T13*T31**2
    MF[6][6]=T11*T33**2 + 2*T13*T31*T33
    MF[6][7]=T12*T33**2 + 2*T13*T32*T33
    MF[6][8]=2*T12*T32*T33 + T13*T32**2
    MF[6][9]=T11*T32*T33 + T12*T31*T33 + T13*T31*T32
    MF[7][0]=3*T21*T31**2
    MF[7][1]=3*T22*T32**2
    MF[7][2]=3*T23*T33**2
    MF[7][3]=T21*T32**2 + 2*T22*T31*T32
    MF[7][4]=2*T21*T31*T32 + T22*T31**2
    MF[7][5]=2*T21*T31*T33 + T23*T31**2
    MF[7][6]=T21*T33**2 + 2*T23*T31*T33
    MF[7][7]=T22*T33**2 + 2*T23*T32*T33
    MF[7][8]=2*T22*T32*T33 + T23*T32**2
    MF[7][9]=T21*T32*T33 + T22*T31*T33 + T23*T31*T32
    MF[8][0]=3*T21**2*T31
    MF[8][1]=3*T22**2*T32
    MF[8][2]=3*T23**2*T33
    MF[8][3]=2*T21*T22*T32 + T22**2*T31
    MF[8][4]=T21**2*T32 + 2*T21*T22*T31
    MF[8][5]=T21**2*T33 + 2*T21*T23*T31
    MF[8][6]=2*T21*T23*T33 + T23**2*T31
    MF[8][7]=2*T22*T23*T33 + T23**2*T32
    MF[8][8]=T22**2*T33 + 2*T22*T23*T32
    MF[8][9]=T21*T22*T33 + T21*T23*T32 + T22*T23*T31
    MF[9][0]=6*T11*T21*T31
    MF[9][1]=6*T12*T22*T32
    MF[9][2]=6*T13*T23*T33
    MF[9][3]=2*T11*T22*T32 + 2*T12*T21*T32 + 2*T12*T22*T31
    MF[9][4]=2*T11*T21*T32 + 2*T11*T22*T31 + 2*T12*T21*T31
    MF[9][5]=2*T11*T21*T33 + 2*T11*T23*T31 + 2*T13*T21*T31
    MF[9][6]=2*T11*T23*T33 + 2*T13*T21*T33 + 2*T13*T23*T31
    MF[9][7]=2*T12*T23*T33 + 2*T13*T22*T33 + 2*T13*T23*T32
    MF[9][8]=2*T12*T22*T33 + 2*T12*T23*T32 + 2*T13*T22*T32
    MF[9][9]=T11*T22*T33 + T11*T23*T32 + T12*T21*T33 + T12*T23*T31 + T13*T21*T32 + T13*T22*T31
    Mr=np.linalg.inv(MF)
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    tcoefs*=np.array(keeps)
    fcoefs=Mi@tcoefs
    return fcoefs

def decomOrbitalG(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    T11=T[0,0]
    T12=T[0,1]
    T13=T[0,2]
    T21=T[1,0]
    T22=T[1,1]
    T23=T[1,2]
    T31=T[2,0]
    T32=T[2,1]
    T33=T[2,2]
    MG=np.zeros((15,15))
    MG[0][0]=T33**4
    MG[0][1]=T32*T33**3
    MG[0][2]=T32**2*T33**2
    MG[0][3]=T32**3*T33
    MG[0][4]=T32**4
    MG[0][5]=T31*T33**3
    MG[0][6]=T31*T32*T33**2
    MG[0][7]=T31*T32**2*T33
    MG[0][8]=T31*T32**3
    MG[0][9]=T31**2*T33**2
    MG[1][0]=4*T23*T33**3
    MG[1][1]=T22*T33**3 + 3*T23*T32*T33**2
    MG[1][2]=2*T22*T32*T33**2 + 2*T23*T32**2*T33
    MG[1][3]=3*T22*T32**2*T33 + T23*T32**3
    MG[1][4]=4*T22*T32**3
    MG[1][5]=T21*T33**3 + 3*T23*T31*T33**2
    MG[1][6]=T21*T32*T33**2 + T22*T31*T33**2 + 2*T23*T31*T32*T33
    MG[1][7]=T21*T32**2*T33 + 2*T22*T31*T32*T33 + T23*T31*T32**2
    MG[1][8]=T21*T32**3 + 3*T22*T31*T32**2
    MG[1][9]=2*T21*T31*T33**2 + 2*T23*T31**2*T33
    MG[2][0]=6*T23**2*T33**2
    MG[2][1]=3*T22*T23*T33**2 + 3*T23**2*T32*T33
    MG[2][2]=T22**2*T33**2 + 4*T22*T23*T32*T33 + T23**2*T32**2
    MG[2][3]=3*T22**2*T32*T33 + 3*T22*T23*T32**2
    MG[2][4]=6*T22**2*T32**2
    MG[2][5]=3*T21*T23*T33**2 + 3*T23**2*T31*T33
    MG[2][6]=T21*T22*T33**2 + 2*T21*T23*T32*T33 + 2*T22*T23*T31*T33 + T23**2*T31*T32
    MG[2][7]=2*T21*T22*T32*T33 + T21*T23*T32**2 + T22**2*T31*T33 + 2*T22*T23*T31*T32
    MG[2][8]=3*T21*T22*T32**2 + 3*T22**2*T31*T32
    MG[2][9]=T21**2*T33**2 + 4*T21*T23*T31*T33 + T23**2*T31**2
    MG[3][0]=4*T23**3*T33
    MG[3][1]=3*T22*T23**2*T33 + T23**3*T32
    MG[3][2]=2*T22**2*T23*T33 + 2*T22*T23**2*T32
    MG[3][3]=T22**3*T33 + 3*T22**2*T23*T32
    MG[3][4]=4*T22**3*T32
    MG[3][5]=3*T21*T23**2*T33 + T23**3*T31
    MG[3][6]=2*T21*T22*T23*T33 + T21*T23**2*T32 + T22*T23**2*T31
    MG[3][7]=T21*T22**2*T33 + 2*T21*T22*T23*T32 + T22**2*T23*T31
    MG[3][8]=3*T21*T22**2*T32 + T22**3*T31
    MG[3][9]=2*T21**2*T23*T33 + 2*T21*T23**2*T31
    MG[4][0]=T23**4
    MG[4][1]=T22*T23**3
    MG[4][2]=T22**2*T23**2
    MG[4][3]=T22**3*T23
    MG[4][4]=T22**4
    MG[4][5]=T21*T23**3
    MG[4][6]=T21*T22*T23**2
    MG[4][7]=T21*T22**2*T23
    MG[4][8]=T21*T22**3
    MG[4][9]=T21**2*T23**2
    MG[5][0]=4*T13*T33**3
    MG[5][1]=T12*T33**3 + 3*T13*T32*T33**2
    MG[5][2]=2*T12*T32*T33**2 + 2*T13*T32**2*T33
    MG[5][3]=3*T12*T32**2*T33 + T13*T32**3
    MG[5][4]=4*T12*T32**3
    MG[5][5]=T11*T33**3 + 3*T13*T31*T33**2
    MG[5][6]=T11*T32*T33**2 + T12*T31*T33**2 + 2*T13*T31*T32*T33
    MG[5][7]=T11*T32**2*T33 + 2*T12*T31*T32*T33 + T13*T31*T32**2
    MG[5][8]=T11*T32**3 + 3*T12*T31*T32**2
    MG[5][9]=2*T11*T31*T33**2 + 2*T13*T31**2*T33
    MG[6][0]=12*T13*T23*T33**2
    MG[6][1]=3*T12*T23*T33**2 + 3*T13*T22*T33**2 + 6*T13*T23*T32*T33
    MG[6][2]=2*T12*T22*T33**2 + 4*T12*T23*T32*T33 + 4*T13*T22*T32*T33 + 2*T13*T23*T32**2
    MG[6][3]=6*T12*T22*T32*T33 + 3*T12*T23*T32**2 + 3*T13*T22*T32**2
    MG[6][4]=12*T12*T22*T32**2
    MG[6][5]=3*T11*T23*T33**2 + 3*T13*T21*T33**2 + 6*T13*T23*T31*T33
    MG[6][6]=T11*T22*T33**2 + 2*T11*T23*T32*T33 + T12*T21*T33**2 + 2*T12*T23*T31*T33 + 2*T13*T21*T32*T33 + 2*T13*T22*T31*T33 + 2*T13*T23*T31*T32
    MG[6][7]=2*T11*T22*T32*T33 + T11*T23*T32**2 + 2*T12*T21*T32*T33 + 2*T12*T22*T31*T33 + 2*T12*T23*T31*T32 + T13*T21*T32**2 + 2*T13*T22*T31*T32
    MG[6][8]=3*T11*T22*T32**2 + 3*T12*T21*T32**2 + 6*T12*T22*T31*T32
    MG[6][9]=2*T11*T21*T33**2 + 4*T11*T23*T31*T33 + 4*T13*T21*T31*T33 + 2*T13*T23*T31**2
    MG[7][0]=12*T13*T23**2*T33
    MG[7][1]=3*T12*T23**2*T33 + 6*T13*T22*T23*T33 + 3*T13*T23**2*T32
    MG[7][2]=4*T12*T22*T23*T33 + 2*T12*T23**2*T32 + 2*T13*T22**2*T33 + 4*T13*T22*T23*T32
    MG[7][3]=3*T12*T22**2*T33 + 6*T12*T22*T23*T32 + 3*T13*T22**2*T32
    MG[7][4]=12*T12*T22**2*T32
    MG[7][5]=3*T11*T23**2*T33 + 6*T13*T21*T23*T33 + 3*T13*T23**2*T31
    MG[7][6]=2*T11*T22*T23*T33 + T11*T23**2*T32 + 2*T12*T21*T23*T33 + T12*T23**2*T31 + 2*T13*T21*T22*T33 + 2*T13*T21*T23*T32 + 2*T13*T22*T23*T31
    MG[7][7]=T11*T22**2*T33 + 2*T11*T22*T23*T32 + 2*T12*T21*T22*T33 + 2*T12*T21*T23*T32 + 2*T12*T22*T23*T31 + 2*T13*T21*T22*T32 + T13*T22**2*T31
    MG[7][8]=3*T11*T22**2*T32 + 6*T12*T21*T22*T32 + 3*T12*T22**2*T31
    MG[7][9]=4*T11*T21*T23*T33 + 2*T11*T23**2*T31 + 2*T13*T21**2*T33 + 4*T13*T21*T23*T31
    MG[8][0]=4*T13*T23**3
    MG[8][1]=T12*T23**3 + 3*T13*T22*T23**2
    MG[8][2]=2*T12*T22*T23**2 + 2*T13*T22**2*T23
    MG[8][3]=3*T12*T22**2*T23 + T13*T22**3
    MG[8][4]=4*T12*T22**3
    MG[8][5]=T11*T23**3 + 3*T13*T21*T23**2
    MG[8][6]=T11*T22*T23**2 + T12*T21*T23**2 + 2*T13*T21*T22*T23
    MG[8][7]=T11*T22**2*T23 + 2*T12*T21*T22*T23 + T13*T21*T22**2
    MG[8][8]=T11*T22**3 + 3*T12*T21*T22**2
    MG[8][9]=2*T11*T21*T23**2 + 2*T13*T21**2*T23
    MG[9][0]=6*T13**2*T33**2
    MG[9][1]=3*T12*T13*T33**2 + 3*T13**2*T32*T33
    MG[9][2]=T12**2*T33**2 + 4*T12*T13*T32*T33 + T13**2*T32**2
    MG[9][3]=3*T12**2*T32*T33 + 3*T12*T13*T32**2
    MG[9][4]=6*T12**2*T32**2
    MG[9][5]=3*T11*T13*T33**2 + 3*T13**2*T31*T33
    MG[9][6]=T11*T12*T33**2 + 2*T11*T13*T32*T33 + 2*T12*T13*T31*T33 + T13**2*T31*T32
    MG[9][7]=2*T11*T12*T32*T33 + T11*T13*T32**2 + T12**2*T31*T33 + 2*T12*T13*T31*T32
    MG[9][8]=3*T11*T12*T32**2 + 3*T12**2*T31*T32
    MG[9][9]=T11**2*T33**2 + 4*T11*T13*T31*T33 + T13**2*T31**2
    Mr=np.linalg.inv(MG)
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    tcoefs*=np.array(keeps)
    fcoefs=Mi@tcoefs
    return fcoefs