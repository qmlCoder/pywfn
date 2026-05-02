from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.geome import Geome
from pywfn.reader import Reader
import re
import numpy as np
from pywfn import core

class WfnReader(Reader):
    def __init__(self, path: str) -> None:
        super().__init__(path)
        match=re.match(r'GAUSSIAN +(\d)+ +MOL ORBITALS +(\d+) PRIMITIVES +(\d+) NUCLEI',self.getline(1))
        assert match is not None, "Invalid WFN file format"
        print(match.groups())
        nobt,nbas,natm=[int(e) for e in match.groups()]
        self.nobt=nobt
        self.nbas=nbas
        self.natm=natm
        self.reader=core.reader.WfnReader(path) # type: ignore

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
        coefs_core=self.reader.get_ciefs()
        coefs=Coefs()
        coefs.core=coefs_core
        return coefs
    
    def read_basAtms(self):
        """读取每个高斯函数对应的原子"""
        basAtms:list[int]=[]
        start=2+self.natm
        count=self.nbas//20 + 1 if self.nbas%20!=0 else 0
        lines=self.getlines(start,start+count)
        for line in lines:
            # print(line)
            atms=line.split()[2:]
            atms=[int(a) for a in atms]
            # print(atms)
            basAtms.extend(atms)
        return basAtms
    
    def read_basSyms(self):
        """读取每个高斯基函数对应的符号"""
        syms=['S','PX','PY','PZ','XX','YY','ZZ','XY','XZ','YZ','XXX','YYY','ZZZ','XXY','XXZ','YYZ','XYY','XZZ','YZZ','XYZ']
        basSyms:list[str]=[]
        count=self.nbas//20 + 1 if self.nbas%20!=0 else 0
        start=2+self.natm+count
        lines=self.getlines(start,start+count)
        for line in lines:
            syms=[syms[int(e)-1] for e in line.split()[2:]]
            basSyms.extend(syms)
        return basSyms
    
    def read_basExps(self):
        """读取每个高斯基函数对应的系数"""
        basExps:list[float]=[]
        count=self.nbas//20 + 1 if self.nbas%20!=0 else 0
        start=2+self.natm+count+count
        lines=self.getlines(start,start+count)
        for line in lines:
            # print(line)
            coefs=[float(e.replace('D','E')) for e in line.split()[1:]]
            # print(coefs)
            basExps.extend(coefs)
        return basExps
    
    def read_coefs(self):
        """读取轨道系数"""
        count0=self.nbas//20 + 1 if self.nbas%20!=0 else 0
        count1=self.nbas//5 + 1 if self.nbas%5!=0 else 0
        start=2+self.natm+count0*2+count1
        Mat=np.zeros(shape=(self.nbas,self.nobt))
        for o in range(self.nobt):
            lines=self.getlines(start,start+count1+1)
            print(start,lines[0])
            vals=[]
            for line in lines[1:]:
                vals+=[float(e.replace('D','E')) for e in line.strip().split()]
            Mat[:,o]=vals
            start+=count1+1
        return Mat