"""
记录基函数变换矩阵
http://sobereva.com/97 感谢卢天大大
"""
import numpy as np
from pywfn import core

DMAT:np.ndarray=core.datas.btrans.DMAT() # type: ignore
FMAT:np.ndarray=core.datas.btrans.FMAT() # type: ignore
GMAT:np.ndarray=core.datas.btrans.GMAT() # type: ignore
HMAT:np.ndarray=core.datas.btrans.HMAT() # type: ignore

MIX_S_LMNS=core.datas.btrans.MIX_S_LMNS() # type: ignore
MIX_P_LMNS=core.datas.btrans.MIX_P_LMNS() # type: ignore
CAR_D_LMNS=core.datas.btrans.CAR_D_LMNS() # type: ignore
CAR_F_LMNS=core.datas.btrans.CAR_F_LMNS() # type: ignore
CAR_G_LMNS=core.datas.btrans.CAR_G_LMNS() # type: ignore
CAR_H_LMNS=core.datas.btrans.CAR_H_LMNS() # type: ignore

Mats:list[np.ndarray]=[None,None,DMAT,FMAT,GMAT,HMAT] # type: ignore

carSsyms=['S']
carPsyms=['PX','PY','PZ']

carDsyms=[
    'XX', # [2,0,0]
    'YY', # [0,2,0]
    'ZZ', # [0,0,2]
    'XY', # [1,1,0]
    'XZ', # [1,0,1]
    'YZ', # [0,1,1]
]
carFsyms=[
    'XXX', # [3,0,0]
    'YYY', # [0,3,0]
    'ZZZ', # [0,0,3]
    'XYY', # [1,2,0]
    'XXY', # [2,1,0]
    'XXZ', # [2,0,1]
    'XZZ', # [1,0,2]
    'YZZ', # [0,1,2]
    'YYZ', # [0,2,1]
    'XYZ', # [1,1,1]
]
carGsyms=[
    'ZZZZ', # [0,0,4]
    'YZZZ', # [0,1,3]
    'YYZZ', # [0,2,2]
    'YYYZ', # [0,3,1]
    'YYYY', # [0,4,0]
    'XZZZ', # [1,0,3]
    'XYZZ', # [1,1,2]
    'XYYZ', # [1,2,1]
    'XYYY', # [1,3,0]
    'XXZZ', # [2,0,2]
    'XXYZ', # [2,1,1]
    'XXYY', # [2,2,0]
    'XXXZ', # [3,0,1]
    'XXXY', # [3,1,0]
    'XXXX', # [3,0,0]
]
carHsyms=[
    'ZZZZZ', # [0,0,5]
    'YZZZZ', # [0,1,4]
    'YYZZZ', # [0,2,3]
    'YYYZZ', # [0,3,2]
    'YYYYZ', # [0,4,1]
    'YYYYY', # [0,5,0]
    'XZZZZ', # [1,0,4]
    'XYZZZ', # [1,1,3]
    'XYYZZ', # [1,2,2]
    'XYYYZ', # [1,3,1]
    'XYYYY', # [1,4,0]
    'XXZZZ', # [2,0,3]
    'XXYZZ', # [2,1,2]
    'XXYYZ', # [2,2,1]
    'XXYYY', # [2,3,0]
    'XXXZZ', # [3,0,2]
    'XXXYZ', # [3,1,1]
    'XXXYY', # [3,2,0]
    'XXXXZ', # [4,0,1]
    'XXXXY', # [4,1,0]
    'XXXXX', # [5,0,0]
]

carSyms=[carSsyms,carPsyms,carDsyms,carFsyms,carGsyms,carHsyms]


sphDsyms=['D 0','D+1','D-1','D+2','D-2']
sphFsyms=['F 0','F+1','F-1','F+2','F-2','F+3','F-3']
sphGsyms=['G 0','G+1','G-1','G+2','G-2','G+3','G-3','G+4','G-4']
sphHsyms=['H 0','H+1','H-1','H+2','H-2','H+3','H-3','H+4','H-4','H+5','H-5']

sphSyms=[carSsyms,carPsyms,sphDsyms,sphFsyms,sphGsyms,sphHsyms]