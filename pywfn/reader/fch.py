"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""
import re
import numpy as np

from pywfn import base
from pywfn.data.elements import elements
from pywfn import reader
from collections import defaultdict
import math

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

        self.jobTitle:str=self.getline(0)
        
        mathc2=re.match(r'(.{10})(.{30})(.{30})',self.getline(1)).groups() # type: ignore
        self.jobType,self.jobMethod,self.jobBasis=mathc2
        self.needs=[
            'Atomic numbers', # 原子种类
            'Current cartesian coordinates', # 原子坐标
            'Shell to atom map',
            'Primitive exponents',
            'Contraction coefficients',
            'Alpha Orbital Energies',
            'Alpha MO coefficients',
            'Total Energy', # 分子能量
            'Charge', # 分子电荷
            'Multiplicity', # 自旋多重度
            'Alpha MO coefficients',
            'Beta MO coefficients',
            'Mulliken Chrgs',
            'Total SCF Density',
            'Orthonormal basis'
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
            
    def parse_title(self,title:str):
        lineNum=self.titles[title]
        line=self.getline(lineNum)
        finds=[line[u:l] for u,l in self.marks]
        _,_,dtype,hasData,dataNum=finds
        dtype=dtype.strip()
        hasData=bool(hasData)
        dataNum=int(dataNum)
        lineSpan=dataNum//6+(0 if dataNum%6==0 else 1)
        tpmap={'I':rf'\d+','R':rf'-?\d+.\d+'}
        text='\n'.join(self.getlines(lineNum+1,lineNum+1+lineSpan+1))
        values=re.findall(tpmap[dtype],text)
        if dtype=='I':
            values=[int(v) for v in values]
        if dtype=='R':
            values=[float(v) for v in values]
        return values
    
    def get_symbols(self) -> list[str]:
        atomics:list[str]=self.read_atoms() # type: ignore
        symbols=[elements[a].symbol for a in atomics]
        return symbols
    
    def get_coords(self) -> np.ndarray:
        coords=self.read_coord()
        return coords
    
    def get_energy(self) -> float:
        lineNum=self.titles['Total Energy']
        line=self.getline(lineNum)
        energy=float(line[45:71])
        return energy
    
    def get_charge(self)->int:
        lineNum=self.titles['Charge']
        line=self.getline(lineNum)
        charge=int(line[45:71])
        return charge
    
    def get_spin(self)->int:
        lineNum=self.titles['Multiplicity']
        line=self.getline(lineNum)
        charge=int(line[45:71])
        return charge
    
    def get_CM(self) -> np.ndarray:
        lineNumA=self.titles['Alpha MO coefficients']
        lineNumB=self.titles['Beta MO coefficients']
        CMs=[]
        for lineNum in [lineNumA,lineNumB]:
            if lineNum==0:continue
            line=self.getline(lineNum)
            # print(line)
            find=re.search(rf'N= +(\d+)',line)
            assert find is not None,"不能找到系数矩阵"
            nval=int(find.groups()[0])
            nmat=int(nval**0.5)
            nline=math.ceil(nval/5)
            lines=self.getlines(lineNumA+1,lineNumA+1+nline)
            vals=re.findall(rf'-?\d+.\d+E[+-]\d+',''.join(lines))
            CM=np.array(vals,dtype=float).reshape(nmat,nmat).T
            CMs.append(CM)
        return np.concatenate(CMs,axis=0)

    def get_OB(self):
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

    def get_DM(self)->np.ndarray:
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
        values=np.array(values,dtype=np.float32).reshape(-1,3)/1.889726
        return values
    
    
