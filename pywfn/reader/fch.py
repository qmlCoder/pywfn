"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""
import re
import numpy as np

from pywfn.base.basis import BasisData,Basis
from pywfn.base.coefs import Coefs
from pywfn.base.geome import Geome
from pywfn.data.elements import elements
from pywfn import reader
from pywfn.reader.utils import toCart
from collections import defaultdict
import math
from functools import lru_cache
from pywfn.data import consts
from pywfn import core


class FchReader(reader.Reader):
    def __init__(self,path:str):
        super().__init__(path)
        self.reader=core.reader.FchReader(path) # type: ignore

    
    def get_geome(self)->"Geome":
        """获取分子几何信息"""
        
        geome_core=self.reader.get_geome()
        geome=Geome()
        geome.core=geome_core
        return geome
    
    def get_basis(self)->"Basis":
        basis_core=self.reader.get_basis()
        basis=Basis()
        basis.core=basis_core
        return basis
    
    
    def get_coefs(self)->"Coefs":
        coefs_core=self.reader.get_coefs()
        coefs=Coefs()
        coefs.core=coefs_core
        return coefs
    
    # def get_energy(self) -> float:
    #     lineNum=self.titles['Total Energy']
    #     line=self.getline(lineNum)
    #     energy=float(line[45:71])
    #     return energy
    
    # def get_nele(self) -> tuple[int, int]:
    #     elea=int(self.getline(self.titles['Number of alpha electrons'])[45:71])
    #     eleb=int(self.getline(self.titles['Number of beta electrons'])[45:71])
    #     return elea,eleb

    # def read_OB(self):
    #     lineNum=self.titles['Orthonormal basis']
    #     line=self.getline(lineNum)
    #     find=re.search(rf'N= +(\d+)',line)
    #     assert find is not None,"不能找到键级矩阵"
    #     nval=int(find.groups()[0])
    #     nmat=int(nval**0.5)
    #     nline=math.ceil(nval/5)
    #     lines=self.getlines(lineNum+1,lineNum+1+nline)
    #     vals=re.findall(rf'-?\d+.\d+E[+-]\d+',''.join(lines))
    #     return np.array(vals,dtype=float).reshape(nmat,nmat).T

    # def read_DM(self)->np.ndarray:
    #     lineNum=self.titles['Total SCF Density']
    #     line=self.getline(lineNum)
    #     find=re.search(rf'N= +(\d+)',line)
    #     assert find is not None,"不能找到密度矩阵"
    #     nval=int(find.groups()[0])
    #     nmat=int(((1+8*nval)**0.5-1)/2)
    #     nline=math.ceil(nval/5)
    #     lines=self.getlines(lineNum+1,lineNum+1+nline)
    #     vals=re.findall(rf'-?\d+.\d+E[+-]\d+',''.join(lines))
    #     DM=np.zeros(shape=(nmat,nmat))
    #     idx=0
    #     for i in range(nmat):
    #         for j in range(i+1):
    #             DM[i,j]=vals[idx]
    #             DM[j,i]=vals[idx]
    #             idx+=1
    #     return DM
    
    # def read_atoms(self):
    #     values=self.parse_title('Atomic numbers')
    #     return values

    # def read_coord(self):
    #     values=self.parse_title('Current cartesian coordinates')
    #     values=np.array(values,dtype=np.float32).reshape(-1,3)
    #     return values
    
    # def read_basis(self):
    #     shlTypes:list[int]=self.parse_title('Shell types') # type: ignore
    #     shlNums:list[int]=self.parse_title('Number of primitives per shell') # type: ignore
    #     shlAtms:list[int]=self.parse_title('Shell to atom map') # type: ignore
    #     pmexps=self.parse_title('Primitive exponents')
    #     pmcoes=self.parse_title('Contraction coefficients')
    #     assert len(pmcoes)!=0,"没有找到系数"
    #     assert len(pmexps)!=0,"没有找到指数"
    #     if self.titles['P(S=P) Contraction coefficients']:
    #         spcoes=self.parse_title('P(S=P) Contraction coefficients')
    #     basisDatas:list[BasisData]=[]
    #     atm=1
    #     shl=1
    #     idx=0
    #     for i in range(len(shlTypes)):
    #         shlType=shlTypes[i]
    #         shlNum=shlNums[i]
    #         shlAtm=shlAtms[i]
            
    #         if atm!=shlAtm:
    #             atm=shlAtm
    #             shl=1
    #         for j in range(shlNum): # type: ignore
    #             if shlType==-1:
    #                 basisDatas.append(BasisData(shlAtm,shl,0,pmexps[idx],pmcoes[idx]))
    #                 basisDatas.append(BasisData(shlAtm,shl,1,pmexps[idx],spcoes[idx]))
    #             else:
    #                 basisDatas.append(BasisData(shlAtm,shl,abs(shlType),pmexps[idx],pmcoes[idx]))
    #             idx+=1
    #         shl+=1
    #     # for each in basisData:
    #     #     print(each)
    #     basisDatas.sort(key=lambda b:(b.atm,b.shl,b.ang))
    #     return basisDatas

    # def read_method(self)->str:
    #     return self.getline(1)[10:20].strip()
    
    # def read_CM(self) -> np.ndarray:
    #     acoefs=self.parse_title('Alpha MO coefficients')
    #     bcoefs=self.parse_title('Beta MO coefficients')
    #     nmat=int(len(acoefs)**0.5)
    #     CMA=np.array(acoefs).reshape(nmat,nmat).T
    #     if len(bcoefs)==0:
    #         return CMA
    #     else:
    #         CMB=np.array(bcoefs).reshape(nmat,nmat).T
    #         return np.concatenate([CMA,CMB],axis=1)
    
    # @lru_cache
    # def read_CMs(self):
    #     shlTypes:list[int]=self.parse_title('Shell types') # type: ignore
    #     shlAtms:list[int]=self.parse_title('Shell to atom map') # type: ignore
    #     atms=[]
    #     shls=[]
    #     syms=[]
    #     shl=1
    #     atm=0
    #     for i,st in enumerate(shlTypes):
    #         assert st in SHL2SYM.keys(),"不认识的shell类型"
    #         syms+=SHL2SYM[st]
    #         atms+=[shlAtms[i]]*len(SHL2SYM[st])
    #         if atm!=shlAtms[i]:
    #             shl=1
    #             atm=shlAtms[i]
    #         else:
    #             shl+=1
    #         shls+=[shl]*len(SHL2SYM[st])
    #     CM=self.read_CM()
    #     return toCart(atms,shls,syms,CM)
