"""
记录基函数变换矩阵
http://sobereva.com/97 感谢卢天大大
"""
import numpy as np

DMat=np.zeros(shape=(6,5))
DMat[ 0, 0]=-0.5
DMat[ 0, 3]=np.sqrt(3)/2
DMat[ 1, 0]=-0.5
DMat[ 1, 3]=-np.sqrt(3)/2    #[1,3]
DMat[ 2, 0]= 1.0
DMat[ 3, 4]= 1.0
DMat[ 4, 1]= 1.0
DMat[ 5, 2]= 1.0

FMat=np.zeros(shape=(10,7))
FMat[ 0, 1]=-np.sqrt(3/8)
FMat[ 0, 5]= np.sqrt(3/8) #5/8
FMat[ 1, 2]=-np.sqrt(3/8)
FMat[ 1, 6]=-np.sqrt(3/8) #5/8
FMat[ 2, 0]= 1.0
FMat[ 3, 1]=-np.sqrt(3/40)
FMat[ 3, 5]=-3/np.sqrt(8) 
FMat[ 4, 2]=-np.sqrt(3/40)
FMat[ 4, 6]=-3/np.sqrt(8)
FMat[ 5, 0]=-3/2/np.sqrt(5)
FMat[ 5, 3]= np.sqrt(3)/2
FMat[ 6, 1]= np.sqrt(6/5)
FMat[ 7, 2]= np.sqrt(6/5)
FMat[ 8, 0]=-3/2/np.sqrt(5)
FMat[ 8, 3]=-np.sqrt(3)/2
FMat[ 9, 4]= 1.0

GMat=np.zeros(shape=(15,9))
GMat[ 0, 0]= 1.0
GMat[ 1, 2]= 2*np.sqrt(5/14)
GMat[ 2, 0]=-3*np.sqrt(3/35)
GMat[ 2, 3]=-3*np.sqrt(3/28)
GMat[ 3, 2]=-3/2*np.sqrt(5/14)
GMat[ 3, 6]=-np.sqrt(5/8)
GMat[ 4, 0]= 3/8
GMat[ 4, 3]= np.sqrt(5)/4
GMat[ 4, 7]= np.sqrt(35)/8
GMat[ 5, 1]= 2*np.sqrt(5/8) #2*np.sqrt(5/14)
GMat[ 6, 4]= 3/np.sqrt(7)
GMat[ 7, 1]=-3/2/np.sqrt(14)
GMat[ 7, 5]=-3/np.sqrt(8)
GMat[ 8, 4]=-np.sqrt(5/28)
GMat[ 8, 8]=-np.sqrt(5/2)
GMat[ 9, 0]=-3*np.sqrt(3/35)
GMat[ 9, 3]= 3*np.sqrt(3/28)
GMat[10, 2]=-3/2*np.sqrt(14)
GMat[10, 6]= 3/np.sqrt(8)
GMat[11, 0]= 3/4*np.sqrt(3/35)
GMat[11, 7]=-3/4*np.sqrt(3)
GMat[12, 1]=-3/2*np.sqrt(5/14)
GMat[12, 5]= np.sqrt(5/8)
GMat[13, 4]=-np.sqrt(5/28)
GMat[13, 8]= np.sqrt(5)/2
GMat[14, 0]= 3/8
GMat[14, 3]=-np.sqrt(5)/4
GMat[14, 7]= np.sqrt(35)/8

HMat=np.zeros(shape=(21,11))
HMat[ 0, 0]= 1.0
HMat[ 1, 2]= np.sqrt(5/3)
HMat[ 2, 0]=-5/np.sqrt(21)
HMat[ 2, 3]=-np.sqrt(5/2)
HMat[ 3, 2]=-3*np.sqrt(5/28)
HMat[ 3, 6]=-np.sqrt(5/6)
HMat[ 4, 0]= 5/8
HMat[ 4, 3]= np.sqrt(35/3)/4
HMat[ 4, 7]= np.sqrt(35/8)
HMat[ 5, 2]= np.sqrt(15/8)
HMat[ 5, 6]= np.sqrt(35/2)/8
HMat[ 5,10]=3/8*np.sqrt(7/2)
HMat[ 6, 1]= np.sqrt(5/3)
HMat[ 7, 4]= np.sqrt(5/3)
HMat[ 8, 1]=-3*np.sqrt(28)
HMat[ 8, 5]=-np.sqrt(3/2)
HMat[ 9, 4]=-np.sqrt(5/12)
HMat[ 9, 8]=-np.sqrt(5/2)
HMat[10, 1]= np.sqrt(5/3)/8
HMat[10, 5]= np.sqrt(35/2)/8
HMat[10, 9]= 5/8*np.sqrt(7/2)
HMat[11, 0]=-5*np.sqrt(21)
HMat[11, 3]= np.sqrt(5/2)
HMat[12, 2]=-3*np.sqrt(28)
HMat[12, 6]= np.sqrt(3/2)
HMat[13, 0]= np.sqrt(15/7)/4
HMat[13, 7]=-3/4*np.sqrt(3)
HMat[14, 2]= np.sqrt(5/7)/4
HMat[14, 6]=-np.sqrt(5/6)/4
HMat[14,10]=-5/4*np.sqrt(3/2)
HMat[15, 1]=-3*np.sqrt(5/28)
HMat[15, 5]= np.sqrt(5/6)
HMat[16, 4]=-np.sqrt(5/12)
HMat[16, 8]= np.sqrt(5/2)
HMat[17, 1]= np.sqrt(5/7)/4
HMat[17, 5]= np.sqrt(5/6)/4
HMat[17, 9]=-5/4*np.sqrt(3/2)
HMat[18, 0]= 5/8
HMat[18, 3]=-np.sqrt(35/3)/4
HMat[18, 7]= np.sqrt(35/8)
HMat[19, 2]= np.sqrt(5/3)/8
HMat[19, 6]=-np.sqrt(35/2)/8
HMat[19,10]=5/8*np.sqrt(7/2)
HMat[20, 1]= np.sqrt(15/8)
HMat[20, 5]=-np.sqrt(35/2)/8
HMat[20, 9]= 3/8*np.sqrt(7/2)

Mats:list[np.ndarray]=[None,None,DMat,FMat,GMat,HMat] # type: ignore

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