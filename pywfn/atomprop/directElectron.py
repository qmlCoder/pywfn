import numpy as np
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures._base import Future

from pywfn.base import Mol,Atom
from pywfn.atomprop import lutils
from pywfn.utils import printer

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.vects:list[np.ndarray]=None
        self.atoms:list[int]=None
        self.zero=True
        self.keep=False
        self.multi=True # 是否开启多线程
        self.abs=True
        
        self.index=0
        self.total=0

    def get_Es(self,obts:list[int],atom:int,vect:np.ndarray,idx:int): 
        """计算电子数量"""
        
        CM_=self.mol.projCM([atom],obts,[vect],zero=self.zero,keep=self.keep,abs=self.abs)
        elect=lutils.get_ects(self.mol,obts,CM_)[atom-1]
        return atom,elect,idx

    def done(self,res:Future):
        atom,elect,idx=res.result()
        self.result[idx]=elect
        self.index+=1
        print(f'\r{self.index}/{self.total}',end='')


    def calculate(self):
        assert self.vects is not None,'未指定方向'
        assert self.atoms is not None,'未指定原子'
        assert len(self.vects)==len(self.atoms),'方向与原子数量不一致'
        self.result=np.zeros(len(self.vects))
        self.total=len(self.vects)
        obts=self.mol.O_obts
        if self.multi: # 多线程计算
            # print('多线程计算')
            self.pool=ThreadPoolExecutor(max_workers=2)
            for i,(atom,vect) in enumerate(zip(self.atoms,self.vects)):
                t=self.pool.submit(self.get_Es,obts,atom,vect,i)
                t.add_done_callback(self.done)
            self.pool.shutdown()
        else:
            # print('单线程计算')
            
            for i,(atom,vect) in enumerate(zip(self.atoms,self.vects)):
                atom,elect,idx=self.get_Es(obts,atom,vect,i)
                self.result[idx]=elect
                # print(f'\r{i}/{self.total}',end='')
        return self.result

    def printRes(self):
        resStr=self.resStr()
        printer.info(f'原子{self.atoms}的方向电子数: ')
        printer.res(resStr)
    
    def resStr(self)->str:
        self.calculate()
        return '\n'.join([f'{a:>4}{e:>10.4f}' for a,e in zip(self.atoms,self.result)])