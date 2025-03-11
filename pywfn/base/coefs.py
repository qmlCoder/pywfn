"""
存储分子轨道系数矩阵
存储直接读到的系数矩阵及相关数据
- CM 系数矩阵
- atoAtms 每一行对应的原子
- atoShls 每一个原子轨道的角动量
- atoSyms 每个原子轨道的符号
- obtEngs 分子轨道能量
- obtOccs 分子轨道是否占据
"""
import numpy as np
from functools import lru_cache,cached_property

from pywfn import base
from pywfn.maths import flib
from pywfn.data import bastrans

def toCart(atms:list[int],shls:list[int],syms:list[str],CM:np.ndarray): # 将数据转为笛卡尔类型的
    CMlist=[]
    atmList=[]
    shlList=[]
    symList=[]
    i=0
    from pywfn.data import bastrans as bt
    while i<CM.shape[0]:
        sym=syms[i]
        match sym:
            case 'D 0':
                CMlist.append(bt.DMat@CM[i:i+5])
                atmList+=[atms[i]]*6
                shlList+=[shls[i]]*6
                symList+=bt.carDsyms
                i+=5
            case 'F 0':
                CMlist.append(bt.FMat@CM[i:i+7])
                atmList+=[atms[i]]*10
                shlList+=[shls[i]]*10
                symList+=bt.carFsyms
                i+=7
            case 'G 0':
                CMlist.append(bt.GMat@CM[i:i+9])
                atmList+=[atms[i]]*15
                shlList+=[shls[i]]*15
                symList+=bt.carGsyms
                i+=9
            case 'H 0':
                CMlist.append(bt.HMat@CM[i:i+11])
                atmList+=[atms[i]]*21
                shlList+=[shls[i]]*21
                symList+=bt.carHsyms
                i+=11
            case 'S'|'PX'|'PY'|'PZ':
                CMlist.append(CM[i][None,:])
                atmList.append(atms[i])
                shlList.append(shls[i])
                symList.append(syms[i])
                i+=1
            case _:
                raise ValueError(f'无效的符号{sym}')

    CM=np.concatenate(CMlist,axis=0)
    return atmList,shlList,symList,CM

class Coefs:
    def __init__(self):
        self.mol:"base.Mol|None"=None
        self._atoAtms_raw:None|list[int]   = None
        self._atoShls_raw:None|list[int]   = None
        self._atoSyms_raw:None|list[str]   = None
        self._obtEngs_raw:None|list[float] = None
        self._obtOccs_raw:None|list[bool]  = None
        self._CM_raw:np.ndarray|None = None
        
        # raw:原始系数，car:笛卡尔系数，sph:球谐系数
    
    @cached_property
    def carData(self)->tuple[list[int],list[int],list[str],np.ndarray]:
        assert self._atoAtms_raw is not None,"未初始化atoAtms"
        assert self._atoShls_raw is not None,"未初始化atoShls"
        assert self._atoSyms_raw is not None,"未初始化atoSyms"
        assert self._CM_raw is not None,"未初始化CM"
        atms,shls,syms,CM=toCart(self._atoAtms_raw,self._atoShls_raw,self._atoSyms_raw,self._CM_raw)
        return atms,shls,syms,CM
    
    @lru_cache
    def CM(self,form:str):
        match form:
            case 'raw':
                assert self._CM_raw is not None,"未初始化CM"
                return self._CM_raw
            case 'car':
                return self.carData[3]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoAtms(self,form:str)->list[int]:
        match form:
            case 'raw':
                assert self._atoAtms_raw is not None,"未初始化atoAtms"
                return self._atoAtms_raw
            case 'car':
                return self.carData[0]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoShls(self,form:str):
        match form:
            case 'raw':
                assert self._atoShls_raw is not None,"未初始化atoShls"
                return self._atoShls_raw
            case 'car':
                return self.carData[1]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoSyms(self,form:str)->list[str]:
        assert self._atoSyms_raw is not None,"未初始化atoSyms"
        match form:
            case 'raw':
                assert self._atoSyms_raw is not None,"未初始化atoSyms"
                return self._atoSyms_raw
            case 'car':
                return self.carData[2]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoKeys(self,form:str):
        atms=self.atoAtms(form)
        shls=self.atoShls(form)
        syms=self.atoSyms(form)
        keys=[]
        for atm,shl,sym in zip(atms,shls,syms):
            keys.append(f'{atm}-{shl}{sym}')
        return keys
    
    def atoXyzs(self,form:str):
        idxs=self.atoAtms(form)
        idxs=[e-1 for e in idxs]
        assert self.mol is not None,"未初始化分子"
        return self.mol.atoms.xyzs[idxs].copy()
    
    @property
    def SM_car(self): # 笛卡尔重叠矩阵
        assert self.mol is not None,"未初始化分子"
        atos,coes,alps,lmns=self.mol.basis.matMap()
        xyzs=self.atoXyzs('car')
        SM=flib.matInteg(atos,coes,alps,lmns,xyzs)
        return SM
    
    @cached_property
    def shlType(self):
        rawTypes={}
        atms=self.atoAtms('raw')
        shls=self.atoShls('raw')
        syms=self.atoSyms('raw') #笛卡尔型符号
        for atm,shl,sym in zip(atms,shls,syms):
            key=f'{atm}-{shl}'
            if key not in rawTypes:rawTypes[key]=''
            if sym in bastrans.carDsyms:
                rawTypes[key]='car'
            elif sym in bastrans.sphDsyms:
                rawTypes[key]='sph'
        return rawTypes
    
    @property
    def SM_raw(self): # 原始重叠矩阵
        # 记录每个原始壳层对应的原始类型
        rawTypes=self.shlType
        SM_car=self.SM_car
        idxs:list[tuple[str,str,int,int]]=[] # 类型，开始索引，结束索引
        atms=self.atoAtms('car')
        shls=self.atoShls('car')
        syms=self.atoSyms('car') #笛卡尔型符号
        idx=0
        # print(syms)
        while True:
            key=f'{atms[idx]}-{shls[idx]}'
            match syms[idx]:
                case 'S':
                    idxs.append((key,'s',idx,idx+1))
                    idx+=1
                case 'PX':
                    idxs.append((key,'p',idx,idx+3))
                    idx+=3
                case 'XX':
                    idxs.append((key,'d',idx,idx+6))
                    idx+=6
                case 'XXX':
                    idxs.append((key,'f',idx,idx+10))
                    idx+=10
                case _:
                    raise ValueError('无效的符号') # 先到这吧
            if idx>=len(syms):break
        # print(idxs)
        SM_mid=[]
        for key,otype,start,end in idxs: #先对行进行替换
            # print(otype,start,end,end-start)
            match otype:
                case 's':
                    SM_mid.append(SM_car[:,start:end])
                case 'p':
                    SM_mid.append(SM_car[:,start:end])
                case 'd':
                    vals=SM_car[:,start:end] # [n,6]
                    if rawTypes[key]=='sph':
                        SM_mid.append(vals@bastrans.DMat)
                    else:
                        SM_mid.append(vals)
                case _:
                    raise ValueError('尚不支持')
        SM_mid=np.concatenate(SM_mid,axis=1)
        SM_raw=[]
        for key,otype,start,end in idxs: #再对列进行替换
            match otype:
                case 's':
                    SM_raw.append(SM_mid[start:end,:])
                case 'p':
                    SM_raw.append(SM_mid[start:end,:])
                case 'd':
                    
                    vals=SM_mid[start:end,:] # [6,n]
                    if rawTypes[key]=='sph':
                        SM_raw.append(bastrans.DMat.T@vals)
                    else:
                        SM_raw.append(vals)
                case _:
                    raise ValueError('尚不支持')
        SM_raw=np.concatenate(SM_raw,axis=0)
        return SM_raw
