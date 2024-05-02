"""
单点能与自由能能量矫正
仅能使用log文件！
"""
from pywfn.reader import LogReader

from pathlib import Path
import numpy as np

class Tool:
    """能量矫正工具"""
    def __init__(self,sfold:str,ffold:str) -> None:
        """
        sfold：single point energy folder 单点能文件所在目录
        ffold：free energy folder 自由能文件所在目录
        """
        self.sfold=sfold
        self.ffold=ffold
        self.sengs=[]
        self.fengs=[]
        self.files=[]

    def extract(self):
        """循环从文件中提取能量"""
        for path in Path(self.sfold).iterdir():
            if path.suffix!='.log':continue
            seng=LogReader(f'{path}').read_energy() # 计算单点能
            self.sengs.append(seng)
            self.files.append(path.name)
        
        for path in Path(self.ffold).iterdir():
            if path.suffix!='.log':continue
            name,engs=LogReader(f'{path}').read_energys()
            feng=engs[3]
            self.fengs.append(feng)

    def save(self,path:str):
        """保存到指定文件"""
        self.extract()
        resStr='file,seng,feng,sum\n'
        for file,seng,feng in zip(self.files,self.sengs,self.fengs):
            resStr+=f'{file},{seng},{feng},{seng+feng}\n'
        Path(path).write_text(resStr)