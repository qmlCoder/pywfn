"""
此脚本用来读取fchk文件
fchk文件中有哪些属性是可以用到的？
"""
import re
import numpy as np

from pywfn import base
from pywfn.data import elements


titleMatch='^(.{40}) {3}(.{1})(.{5})(.{12})$'
class FchReader:
    def __init__(self,path:str):
        self.mol=base.Mol()
        self.path=path
        self.contents:list[Content]=[]
        with open(self.path,'r',encoding='utf-8') as f:
            self.text=f.read()
            self.lines=self.text.split('\n')
        self.jobTitle:str=self.lines[0]
        self.jobType,self.jobMethod,self.jobBasis=re.match(r'(.{10})(.{30})(.{30})',self.lines[1]).groups()
        # 正式的数据从第十行开始
        for idx in range(10,len(self.lines)):
            line=self.lines[idx]
            if re.match(titleMatch, line) is not None:
                self.contents.append(Content(idx+1, line))
        
        self.get_all()
    
    def get_content(self,title:str):
        title=title.ljust(40,' ')
        for each in self.contents:
            if each.title==title:
                return each.get(self)
        else:
            return None
    
    def get_all(self):
        """
        获取所有数据
        1. 原子符号
        2. 原子坐标(fch中的原子坐标单位是Bohr,1Bohr=0.53A)
        3. 轨道能量(根据是否含有alpha判断是否为开壳层)
        4. 系数矩阵
        5. 重叠矩阵
        """
        atoms=self.get_content('Atomic numbers')
        symbols=[elements[int(e)].symbol for e in atoms]
        coords=self.get_content('Current cartesian coordinates').reshape(-1,3)*0.5291772083
        ACMs=self.get_content('Alpha MO coefficients')
        aeN=self.get_content('Number of alpha electrons') #alpha和beta电子数量
        beN=self.get_content('Number of beta electrons')
        basiN=self.get_content('Number of basis functions') # 基函数数量

        for i,atomic in enumerate(atoms):
            symbol=symbols[i]
            coord=list(coords[i])
            self.mol.add_atom(symbol=symbol,coord=coord)
        self.mol.obtElcts=[1]*aeN+[0]*(basiN-aeN)+[1]*beN+[0]*(basiN-beN)
        



class Content:
    def __init__(self,idx:int,line:str):
        """fch中的数据格式非常规整，刚初始化的时候不读取数据，需要的时候再读取"""
        self.idx=idx
        title,dataType,isArray,number=re.match(titleMatch,line).groups()
        self.title:str=title
        self.dataType:str=dataType
        self.isArray=True if 'N=' in isArray else False
        self.number:int=int(number)

    def get(self,reader:FchReader):
        """
        读取其中内容
        整数一行六个数据，浮点数一行五个数据
        """
        ps={
            'I':'\d+',
            'R':'-?\d+.\d+E[\+-]\d{2}'
        }
        ds={
            'I':np.int8,
            'R':np.float32
        }
        lineCount=5 if self.dataType=='R' else 6 #每一行的数据量
        if self.isArray:
            lineNum=self.number//lineCount+1
            text='\n'.join(reader.lines[i] for i in range(self.idx,self.idx+lineNum))
            res=re.findall(ps[self.dataType], text)
            return np.array(res,dtype=ds[self.dataType])

    def __repr__(self):
        return f'{self.idx},{self.title},{self.dataType},{self.isArray},{self.number}\n'
