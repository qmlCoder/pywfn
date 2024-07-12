"""
该脚本用来生成SI(支持信息)文件,包含
能量
坐标
频率
"""
from pathlib import Path
import re
import os

from pywfn.base import Mol
from pywfn.utils import printer
from pywfn.data import temps
from pywfn.reader import LogReader

class Tool:
    def __init__(self,paths:list[str]) -> None:
        """
        读取文件貌似只需要log/out文件就行了
        """
        self.paths=paths
        self.template=temps.si
        self.selects:list[int]=[1,2,3]
        self.sameFile:bool=True
        self.FILENAME:str=''
    
    def write_energy(self,temp:str,reader:LogReader):
        engNums=reader.read_energys()
        """匹配并保存各种校正能量"""
        for i,eng in enumerate(engNums):
            temp=temp.replace(f'<E{i}>',f'{eng:>15.8f}')
        return temp
            
    def write_coords(self,temp:str,reader:LogReader):
        COORDS=[]
        syms,xyzs=reader.read_coords()
        for sym,(x,y,z) in zip(syms,xyzs):
            COORDS.append(f'{sym:<12}{x:16.8f}{y:16.8f}{z:16.8f}')
        temp=temp.replace('<COORD>','\n'.join(COORDS))
        return temp

    def write_freqs(self,temp:str,reader:LogReader):
        freqs=reader.read_freqs()
        freqStr=''
        for i,freq in enumerate(freqs):
            freqStr+=f'{freq:>20.2f}'
            if (i+1)%3==0 and (i+1)!=len(freqs):freqStr+='\n'
        temp=temp.replace('<FREQ>',freqStr)
        return temp
    
    def write_names(self,temp:str,reader:LogReader):
        temp=temp.replace('<NAME>',reader.fname)
        return temp
        # reader.fname

    def build(self):
        """
        导出分子的各种信息到txt文件
        1. 坐标 coord
        2. 能量 energy
        3. 频率 freq
        """
        texts=[]
        for path in self.paths:
            reader=LogReader(path)
            temp=temps.si
            temp=self.write_coords(temp,reader)
            temp=self.write_energy(temp,reader)
            temp=self.write_freqs(temp,reader)
            temp=self.write_names(temp,reader)
            texts.append(temp)
        return texts
        
    def save(self,path:str):
        texts=self.build()
        marke='-'*60+'\n\n' # 分割线
        if Path(path).exists():os.remove(path)
        with open(path,'a',encoding='utf-8') as f:
            f.write(marke.join(texts))
        printer.info(f"SI成功导出 {path} >_<\n")