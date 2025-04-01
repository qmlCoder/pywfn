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
            case _:
                CMlist.append(CM[i][None,:])
                atmList.append(atms[i])
                shlList.append(shls[i])
                symList.append(syms[i])
                i+=1

    CM=np.concatenate(CMlist,axis=0)
    return atmList,shlList,symList,CM

class Coefs:
    def __init__(self):
        self.mol:"base.Mole|None"=None
        self._atoAtms:None|list[int]   = None
        self._atoShls:None|list[int]   = None
        self._atoSyms:None|list[str]   = None
        self.obtEngs:None|list[float] = None
        self.obtOccs:None|list[bool]  = None # 根据α和β电子数计算得到
        self._CM:np.ndarray|None = None
        
        # raw:原始系数，car:笛卡尔系数，sph:球谐系数
    
    @cached_property
    def carData(self)->tuple[list[int],list[int],list[str],np.ndarray]:
        assert self._atoAtms is not None,"未初始化atoAtms"
        assert self._atoShls is not None,"未初始化atoShls"
        assert self._atoSyms is not None,"未初始化atoSyms"
        assert self._CM is not None,"未初始化CM"
        atms,shls,syms,CM=toCart(self._atoAtms,self._atoShls,self._atoSyms,self._CM)
        return atms,shls,syms,CM
    
    @lru_cache
    def CM(self,form:str):
        match form:
            case 'raw':
                assert self._CM is not None,"未初始化CM"
                return self._CM
            case 'car':
                return self.carData[3]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoAtms(self,form:str)->list[int]:
        match form:
            case 'raw':
                assert self._atoAtms is not None,"未初始化atoAtms"
                return self._atoAtms
            case 'car':
                return self.carData[0]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoShls(self,form:str):
        match form:
            case 'raw':
                assert self._atoShls is not None,"未初始化atoShls"
                return self._atoShls
            case 'car':
                return self.carData[1]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    def atoSyms(self,form:str)->list[str]:
        assert self._atoSyms is not None,"未初始化atoSyms"
        match form:
            case 'raw':
                assert self._atoSyms is not None,"未初始化atoSyms"
                return self._atoSyms
            case 'car':
                return self.carData[2]
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')
    
    # @property
    # def obtOccs(self)->list[bool]:
    #     assert self.mol is not None,"未初始化分子"
    #     nela,nelb=self.mol.nele
    #     if self.mol.open:
    #         nobt=self.mol.CM.shape[1]//2
    #         return [True]*nela+[False]*(nobt-nela)+[True]*nelb+[False]*(nobt-nelb)
    #     else:
    #         nobt=self.mol.CM.shape[1]
    #         return [True]*nela+[False]*(nobt-nela)
    
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
    def SM_car_(self): # 笛卡尔重叠矩阵
        from pywfn.maths import flib
        assert self.mol is not None,"未初始化分子"
        atos,coes,alps,lmns=self.mol.basis.basMap()
        xyzs=self.atoXyzs('car')
        facs = [1., 1., 3.]
        nbas=len(coes)
        for i in range(nbas):
            l,m,n=lmns[i]
            fac = facs[l]*facs[m]*facs[n]
            ang=l+m+n
            alp=alps[i]
            Nm=(2.*alp/np.pi)**0.75*np.sqrt((4.*alp)**ang/fac)
            coes[i]=Nm*coes[i]
        SM=flib.matInteg(atos,coes,alps,lmns,xyzs)
        return SM
    
    @property
    def SM_car(self):
        from pywfn.maths import rlib
        assert self.mol is not None,"未初始化分子"
        atos,coes,alps,lmns=self.mol.basis.basMap()
        xyzs=self.atoXyzs('car')
        facs = [1., 1., 3.]
        for i in range(len(coes)):
            l,m,n=lmns[i]
            fac = facs[l]*facs[m]*facs[n]
            ang=l+m+n
            alp=alps[i]
            Nm=(2.*alp/np.pi)**0.75*np.sqrt((4.*alp)**ang/fac)
            coes[i]=Nm*coes[i]
        SM=rlib.mat_integ_rs(atos,coes,alps,lmns,xyzs) # type: ignore
        SM=np.array(SM)
        return SM

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

    def SM(self,form:str):
        match form:
            case 'raw':
                return self.SM_raw
            case 'car':
                return self.SM_car
            case 'sph':
                raise NotImplementedError('球谐系数尚未实现')
            case _:
                raise ValueError('无效的系数形式')

    def CMr(self,form:str): # 行平方和为1的系数矩阵
        SM=self.SM(form)
        CM=self.CM(form)
        eigVal,eigVec=np.linalg.eigh(SM)
        X=(eigVec@np.diag(eigVal))@eigVec.T
        return X.T@CM
    
    def __repr__(self) -> str:
        text=''
        nato,nobt=self.CM('raw').shape
        assert self.obtEngs is not None,"未初始化obtEngs"
        assert self.obtOccs is not None,"未初始化obtOccs"
        assert self._atoAtms is not None,"未初始化_atoAtms"
        assert self._atoShls is not None,"未初始化_atoShls"
        assert self._atoSyms is not None,"未初始化_atoSyms"
        engs=self.obtEngs
        occs=['O' if e else 'V' for e in self.obtOccs]
        atms=self._atoAtms
        shls=self._atoShls
        syms=self._atoSyms
        
        obtIdxs=[0,1,2,3,4]+[nobt-(5-e) for e in range(5)]
        atoIdxs=[0,1,2,3,4]+[nato-(5-e) for e in range(5)]
        text+=' '*15
        for i,obt in enumerate(obtIdxs):
            text+=f'{occs[obt]:>10}'
            if i==4:text+=f"{'.....':>10}"
        text+='\n'
        text+=' '*15
        for i,obt in enumerate(obtIdxs):
            text+=f'{engs[obt]:>10.4f}'
            if i==4:text+=f"{'.....':>10}"
        text+='\n'
        CM=self.CM('raw')
        for j,ato in enumerate(atoIdxs):
            if j==4:text+=(" "*15+"     ....."*11+'\n')
            text+=f"{atms[ato]:>5}{shls[ato]:>5}{syms[ato]:>5}"
            for i,obt in enumerate(obtIdxs):
                text+=f'{CM[ato,obt]:>10.4f}'
                if i==4:text+=f"{'.....':>10}"
            text+='\n'

        return text