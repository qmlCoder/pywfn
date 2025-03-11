"""
将分子对象保存为gif文件
"""
from pathlib import Path
import numpy as np
from pywfn.base import Mol
from pywfn import config
from pywfn.data import temps

class GjfWriter:
    def __init__(self):
        self.syms:list[str]=[]    # 元素符号
        self.xyzs:np.ndarray|None=None # 原子坐标
        self.temp=temps.get('gjf')# 模板
        self.chk='chk'
        self.charge=0
        self.spin  =1
        self.title=config.GJF_TITLE
    
    def fromMol(self,mol:Mol):
        self.syms=mol.atoms.syms
        self.xyzs=mol.atoms.xyzs
        return self
        
    def get_coordStr(self):
        """生成坐标的字符串形式"""
        assert self.xyzs is not None,'请先设置坐标'
        coordStrs=[]
        for i,sym in enumerate(self.syms):
            x,y,z=self.xyzs[i]/config.BOHR_RADIUS
            coordStr=f'{sym:>2}{x:>14.8f}{y:>14.8f}{z:>14.8f}'
            coordStrs.append(coordStr)
        return '\n'.join(coordStrs)
    
    def build(self)->str:
        assert self.title!='','请设置标题'
        self.temp=self.temp.replace('<COORD>',self.get_coordStr())
        self.temp=self.temp.replace('<CHARGE>',f'{self.charge}')
        self.temp=self.temp.replace('<SPIN>',f'{self.spin}')
        self.temp=self.temp.replace('<CHK>',self.chk)
        self.temp=self.temp.replace('<TITLE>',self.title)
        return self.temp

    def save(self,path:str):
        self.build()
        Path(path).write_text(self.temp)
        print(f'文件导出至{path}')