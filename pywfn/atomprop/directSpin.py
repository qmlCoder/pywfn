"""
此脚本用来计算pi-Mulliken电子自旋
"""
import numpy as np
from threading import Thread
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures._base import Future
from pywfn.base import Mol,Atom
from pywfn.utils import printer
from pywfn.atomprop import lutils

class Calculator:
    def __init__(self,mol:"Mol") -> None:
        self.mol=mol
        self.vects:list[np.ndarray]=None
        self.atoms:list[int]=None
        self.zero=True
        self.keep=True
        self.abs=False
        self.index=0
        self.total=0
    
    def get_Es(self,obts:list[int],atom:int,vect:np.ndarray,idx:int,spin:int): 
        """计算电子数量"""
        
        CM_=self.mol.projCM([atom],obts,[vect],zero=self.zero,keep=self.keep,abs=self.abs)
        elect=lutils.get_ects(self.mol,obts,CM_)[atom-1]
        return atom,elect,idx,spin
    
    def cal_group(self,obts:list[int],spin:int):
        for i,(atom,vect) in enumerate(zip(self.atoms,self.vects)):
            t=self.pool.submit(self.get_Es,obts,atom,vect,i,spin)
            t.add_done_callback(self.done)
        
    
    def done(self,res:Future):
        atom,elect,idx,spin=res.result()
        self.result[spin,idx]=elect
        self.index+=1
        print(f'\r{self.index}/{self.total}',end='')

    def calculate(self):
        """计算所有原子的自旋"""
        if not self.mol.isOpenShell:
            printer.warn('非开壳层分子无法计算自旋')
            return 0
        # 首先要有alpha电子和beta电子对应的轨道系数
        self.result=np.zeros(shape=(2,len(self.vects)))
        self.total=len(self.vects)*2
        obtNum=self.mol.CM.shape[0] # 系数矩阵行数，基函数数量
        a_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i< obtNum)] # alpha 电子所在的轨道
        b_obt=[i for i,e in enumerate(self.mol.obtEcts) if (e!=0 and i>=obtNum)] # beta  电子所在的轨道
        # aEs=np.array(self.get_Es(a_obt)) # alpha 电子数量
        # bEs=np.array(self.get_Es(b_obt)) # beta  电子数量
        # return (aEs-bEs)[self.atom-1]
        self.pool=ThreadPoolExecutor(max_workers=10)
        self.cal_group(a_obt,0)
        self.cal_group(b_obt,1)
        self.pool.shutdown()
        # printer.console.log(self.result)
        return self.result[0]-self.result[1]
        
    def printRes(self):
        resStr=self.resStr()
        printer.info('方向电子自旋分布: ')
        printer.res(resStr)
        
    def resStr(self):
        elect=self.calculate()
        return lutils.atomValueStr(self.mol,[self.atom],[elect])