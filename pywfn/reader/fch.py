"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""
import re
import numpy as np

from pywfn import base
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

SHL2SYM={
    -2:['D 0','D+1','D-1','D+2','D-2'],
    -1:['S','PX','PY','PZ'],
    0: ['S'],
    1: ['PX','PY','PZ'],
    2: ['XX','YY','ZZ','XY','XZ','YZ']
}

class Title:
    def __init__(self,lineNum:int,dataNum:int,dataType:str,hasData:bool) -> None:
        self.lineNum=lineNum
        self.dataNum=dataNum
        self.dType=dataType #数据类型
        self.hasData=hasData # 是否包含更多数据

titleMatch='^(.{40}) {3}(.{1})(.{5})(.{12})$'
class FchReader(reader.Reader):
    def __init__(self,path:str):
        super().__init__(path)
        self.type='fch'
        self.jobTitle:str=self.getline(0)
        
        mathc2=re.match(r'(.{10})(.{30})(.{30})',self.getline(1)).groups() # type: ignore
        self.jobType,self.jobMethod,self.jobBasis=mathc2
        self.needs=[
            'Atomic numbers', # 原子种类
            'Current cartesian coordinates', # 原子坐标
            'Primitive exponents',
            'Contraction coefficients',
            'Alpha Orbital Energies',
            'Total Energy',              # 分子能量
            'Alpha MO coefficients',
            'Beta MO coefficients',
            'Number of alpha electrons', # alpha电子数
            'Number of beta electrons', # beta电子数
            'Mulliken Chrgs',
            'Total SCF Density',
            'Orthonormal basis',
            'Shell types', # 壳层类型
            'Number of primitives per shell',
            'Shell to atom map',
            'Primitive exponents', # 基函数指数
            'Contraction coefficients',
            'P(S=P) Contraction coefficients',
            'Alpha Orbital Energies',
            'Beta Orbital Energies'

        ]
        self.titles:dict[str,int]={k:0 for k in self.needs}
        self.marks=[
            [0,40],
            [40,43],
            [43,45],
            [45,50],
            [50,62]
        ]
        self.search_title()

    def search_title(self):
        """搜索需要的数据的标题所在的行"""
        patern=r'[A-Za-z -/]{40}.{3}[IRc]{1}.{5}[ \d]{12}'
        lineNum=self.lineNum
        for l in range(lineNum):
            line=self.getline(l)
            if line[0]==' ':continue
            title=line[0:40].strip()
            if title not in self.needs:continue
            self.titles[title]=l
    
    @lru_cache
    def parse_title(self,title:str)->list[int]|list[float]: # 解析区块内容
        lineNum=self.titles[title]
        if lineNum==0:return []
        line=self.getline(lineNum)
        finds=[line[u:l] for u,l in self.marks]
        _,_,dtype,hasData,dataNum=finds # 数据类型，是否包含更多数据，数据数量
        dtype=dtype.strip()
        hasData=bool(hasData)
        dataNum=int(dataNum)  # 数据的量
        if dtype=='I':
            lineSpan=dataNum//6+(0 if dataNum%6==0 else 1)
            text='\n'.join(self.getlines(lineNum+1,lineNum+1+lineSpan))
            vals=re.findall(rf'-?\d+',text)
            vals=[int(v) for v in vals]
            # print(text)
        if dtype=='R':
            lineSpan=dataNum//5+(0 if dataNum%5==0 else 1)
            text='\n'.join(self.getlines(lineNum+1,lineNum+1+lineSpan))
            vals=re.findall(rf'-?\d+\.\d+E[+-]\d+',text)
            vals=[float(v) for v in vals]
            # print(text)
        return vals
    
    def get_geome(self) -> Geome:
        atomics:list[str]=self.read_atoms() # type: ignore
        syms=[elements[a].symbol for a in atomics]
        xyzs=self.read_coord()
        return Geome().build(syms,xyzs)
    
    def get_energy(self) -> float:
        lineNum=self.titles['Total Energy']
        line=self.getline(lineNum)
        energy=float(line[45:71])
        return energy
    
    def get_nele(self) -> tuple[int, int]:
        elea=int(self.getline(self.titles['Number of alpha electrons'])[45:71])
        eleb=int(self.getline(self.titles['Number of beta electrons'])[45:71])
        return elea,eleb
    
    def get_basis(self)->Basis:
        basData=self.read_basis()
        basName=self.getline(1)[40:90].strip()
        basis=Basis()
        basis.name=basName
        basis.data=basData
        return basis
    
    def get_coefs(self) -> Coefs:
        atms,shls,syms,CM=self.read_CMs()
        aengs=self.parse_title('Alpha Orbital Energies')
        bengs=self.parse_title('Beta Orbital Energies')
        engs:list[float] = aengs+bengs # type: ignore
        coefs=Coefs()
        coefs._atoAtms=atms
        coefs._atoShls=shls
        coefs._atoSyms=syms
        coefs._CM=CM
        coefs.obtEngs=engs
        elea,eleb=self.get_nele()
        coefs.obtOccs=[True]*elea+[False]*(len(aengs)-elea)+[True]*eleb+[False]*(len(bengs)-eleb)
        return coefs

    def read_OB(self):
        lineNum=self.titles['Orthonormal basis']
        line=self.getline(lineNum)
        find=re.search(rf'N= +(\d+)',line)
        assert find is not None,"不能找到键级矩阵"
        nval=int(find.groups()[0])
        nmat=int(nval**0.5)
        nline=math.ceil(nval/5)
        lines=self.getlines(lineNum+1,lineNum+1+nline)
        vals=re.findall(rf'-?\d+.\d+E[+-]\d+',''.join(lines))
        return np.array(vals,dtype=float).reshape(nmat,nmat).T

    def read_DM(self)->np.ndarray:
        lineNum=self.titles['Total SCF Density']
        line=self.getline(lineNum)
        find=re.search(rf'N= +(\d+)',line)
        assert find is not None,"不能找到密度矩阵"
        nval=int(find.groups()[0])
        nmat=int(((1+8*nval)**0.5-1)/2)
        nline=math.ceil(nval/5)
        lines=self.getlines(lineNum+1,lineNum+1+nline)
        vals=re.findall(rf'-?\d+.\d+E[+-]\d+',''.join(lines))
        DM=np.zeros(shape=(nmat,nmat))
        idx=0
        for i in range(nmat):
            for j in range(i+1):
                DM[i,j]=vals[idx]
                DM[j,i]=vals[idx]
                idx+=1
        return DM
    
    def read_atoms(self):
        values=self.parse_title('Atomic numbers')
        return values

    def read_coord(self):
        values=self.parse_title('Current cartesian coordinates')
        values=np.array(values,dtype=np.float32).reshape(-1,3)
        return values
    
    def read_basis(self):
        shlTypes:list[int]=self.parse_title('Shell types') # type: ignore
        shlNums:list[int]=self.parse_title('Number of primitives per shell') # type: ignore
        shlAtms:list[int]=self.parse_title('Shell to atom map') # type: ignore
        pmexps=self.parse_title('Primitive exponents')
        pmcoes=self.parse_title('Contraction coefficients')
        assert len(pmcoes)!=0,"没有找到系数"
        assert len(pmexps)!=0,"没有找到指数"
        if self.titles['P(S=P) Contraction coefficients']:
            spcoes=self.parse_title('P(S=P) Contraction coefficients')
        basisDatas:list[BasisData]=[]
        atm=1
        shl=1
        idx=0
        for i in range(len(shlTypes)):
            shlType=shlTypes[i]
            shlNum=shlNums[i]
            shlAtm=shlAtms[i]
            
            if atm!=shlAtm:
                atm=shlAtm
                shl=1
            for j in range(shlNum): # type: ignore
                if shlType==-1:
                    basisDatas.append(BasisData(shlAtm,shl,0,pmcoes[idx],pmexps[idx]))
                    basisDatas.append(BasisData(shlAtm,shl,1,spcoes[idx],pmexps[idx]))
                else:
                    basisDatas.append(BasisData(shlAtm,shl,abs(shlType),pmcoes[idx],pmexps[idx]))
                idx+=1
            shl+=1
        # for each in basisData:
        #     print(each)
        basisDatas.sort(key=lambda b:(b.atm,b.shl,b.ang))
        return basisDatas

    def read_method(self)->str:
        return self.getline(1)[10:20].strip()
    
    def read_CM(self) -> np.ndarray:
        acoefs=self.parse_title('Alpha MO coefficients')
        bcoefs=self.parse_title('Beta MO coefficients')
        nmat=int(len(acoefs)**0.5)
        CMA=np.array(acoefs).reshape(nmat,nmat).T
        if len(bcoefs)==0:
            return CMA
        else:
            CMB=np.array(bcoefs).reshape(nmat,nmat).T
            return np.concatenate([CMA,CMB],axis=1)
    
    @lru_cache
    def read_CMs(self):
        shlTypes:list[int]=self.parse_title('Shell types') # type: ignore
        shlAtms:list[int]=self.parse_title('Shell to atom map') # type: ignore
        atms=[]
        shls=[]
        syms=[]
        shl=1
        atm=0
        for i,st in enumerate(shlTypes):
            assert st in SHL2SYM.keys(),"不认识的shell类型"
            syms+=SHL2SYM[st]
            atms+=[shlAtms[i]]*len(SHL2SYM[st])
            if atm!=shlAtms[i]:
                shl=1
                atm=shlAtms[i]
            else:
                shl+=1
            shls+=[shl]*len(SHL2SYM[st])
        CM=self.read_CM()
        return toCart(atms,shls,syms,CM)
