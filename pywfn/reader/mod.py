"""
读取.molden文件的脚本
"""
from pywfn import reader
from pywfn.base.basis import BasisData,Basis
from pywfn.base.coefs import Coefs
from pywfn.reader.utils import toCart
from pywfn.data import bastrans
import numpy as np
import re
from functools import lru_cache


sym2ang={'s':0,'p':1,'d':2,'f':3,'g':4,'h':5,'i':6}

class ModReader(reader.Reader):
    def __init__(self, path:str):
        super().__init__(path)
        self.type='mod'
        self.titles={
            'N_ATOMS':0,
            'ATOMS':0,
            'GTO':0,
            'MO':0
        }
        self.basType=''
        self.search_title()
    
    def search_title(self):
        for l in range(self.lineNum):
            line=self.getline(l,keepEnd=False)
            # if line[0]!='[':continue
            if line=='[N_ATOMS]':
                self.titles['N_ATOMS']=l
            elif line=='[Atoms]  AU': # 有一定的风险，不知道molden的单位是不是统一的
                self.titles['ATOMS']=l
            elif line=='[GTO]':
                self.titles['GTO']=l
            elif line=='[MO]':
                self.titles['MO']=l
            elif line=='[5D7F9G]':
                self.basType='sph'
        assert self.basType!='',"不确定基组类型"
    
    def get_nele(self):
        engs,occs,CM=self.read_coefs()
        row,col=CM.shape
        if row>=col: # 闭壳层
            nela=nelb=sum(occs)
        else:
            nela=sum(occs[:col//2])
            nleb=sum(occs[col//2:])
        return nela,nelb
    
    

    def get_atmSyms(self) -> list[str]:
        syms,xyzs=self.read_geom()
        return syms
    
    def get_atmXyzs(self) -> np.ndarray:
        syms,xyzs=self.read_geom()
        return xyzs
    
    def get_basis(self) -> Basis:
        atms,shls,syms,datas=self.read_basis()
        basis=Basis()
        basis.data=datas
        basis.name='mod'
        return basis
    
    def get_coefs(self)->Coefs: # 获取轨道信息
        coefs=Coefs()
        engs,occs,CM=self.read_coefs()
        atms,shls,syms,_=self.read_basis()
        coefs._atoAtms=atms
        coefs._atoShls=shls
        coefs._atoSyms=syms
        coefs._obtEngs=engs
        coefs._CM=CM
        return coefs

    def read_geom(self):
        natm=int(self.getline(self.titles['N_ATOMS']+1).strip()) #原子的数量
        iatm=self.titles['ATOMS']+1
        syms=[]
        xyzs=[]
        for l in range(iatm,iatm+natm):
            line=self.getline(l,keepEnd=False).split()
            # print(line)
            sym,_,_,x,y,z=line
            syms.append(sym)
            xyzs.append([float(x),float(y),float(z)])
        return syms,np.array(xyzs)

    @lru_cache
    def read_basis(self):
        pat1=rf'^ +(\d+)$'
        pat2=rf' +([spdf]) +(\d+)'
        pat3=rf'^ +(-?\d+.\d+E[+-]\d+) +(-?\d+.\d+E[+-]\d+)'
        atm=0
        shl=0
        sym=''
        basisDatas=[]
        syms=['s','p','d','f']
        
        
        atoAtms=[]
        atoShls=[]
        atoSyms=[]
        for l in range(self.titles['GTO']+1,self.titles['MO']-1):
            line=self.getline(l,keepEnd=False)
            if line=='':continue
            if m1:=re.match(pat1,line):
                # print(m1.groups())
                atm=int(m1.groups()[0])
                shl=0
            elif m2:=re.match(pat2,line):
                # print(m2.groups())
                sym,_=m2.groups()
                shl+=1
                match sym,self.basType:
                    case 's',_:
                        atoSyms+=['S']
                        atoAtms+=[atm]*1
                        atoShls+=[shl]*1
                    case 'p',_:
                        atoSyms+=['PX','PY','PZ']
                        atoAtms+=[atm]*3
                        atoShls+=[shl]*3
                    case 'd','sph':
                        atoSyms+=bastrans.sphDsyms
                        atoAtms+=[atm]*5
                        atoShls+=[shl]*5
                    case 'd','car':
                        atoSyms+=bastrans.carDsyms
                        atoAtms+=[atm]*6
                        atoShls+=[shl]*6
                    case 'f','sph':
                        atoSyms+=bastrans.sphFsyms
                        atoAtms+=[atm]*7
                        atoShls+=[shl]*7
                    case 'f','car':
                        atoSyms+=bastrans.carFsyms
                        atoAtms+=[atm]*10
                        atoShls+=[shl]*10
                    
                
            elif m3:=re.match(pat3,line):
                # print(m3.groups())
                alp,coe=m3.groups()
                alp=float(alp)
                coe=float(coe)
                basisDatas.append(BasisData(atm,shl,syms.index(sym),coe,alp))
                
            else:
                print('解析失败',line)
        return atoAtms,atoShls,atoSyms,basisDatas
    
    @lru_cache
    def read_coefs(self): # 读取轨道系数
        pat2=r' Ene= +(-?\d+.\d+)'
        pat4=r' Occup= +(\d.\d+)'
        pat5=r'^ *(\d+) +(-?\d+.\d+)$'
        obtOccs=[]
        obtEngs=[]
        # coes=[]
        CM:list[list[float]]=[]
        for l in range(self.titles['MO']+1,self.lineNum):
            line=self.getline(l,keepEnd=False)
            if line=='':continue
            if line[:5]==' Sym=':continue
            if line[:6]==' Spin=':continue
            if m2:=re.match(pat2,line):
                # print(m2.groups())
                eng=float(m2.groups()[0])
                obtEngs.append(eng)
                CM.append([])
            elif m4:=re.match(pat4,line):
                # print(m4.groups())
                occ=float(m4.groups()[0])
                obtOccs.append(occ!=0.0)
            elif m5:=re.match(pat5,line):
                # print(m5.groups())
                idx,coe=m5.groups()
                CM[-1].append(float(coe))
            else:
                print('解析失败',line)
        return obtEngs,obtOccs,np.array(CM).T