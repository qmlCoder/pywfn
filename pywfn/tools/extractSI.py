"""
该脚本用来生成SI(支持信息)文件,包含
能量
坐标
频率
"""
from pathlib import Path
import re

from pywfn.base import Mol
from pywfn.utils import printer
from pywfn.data import temps
from pywfn.reader import LogReader

class Tool:
    def __init__(self,mol:Mol) -> None:
        """
        读取文件貌似只需要log/out文件就行了
        """
        self.mol=mol
        self.reader:LogReader=mol.reader #已经实例化的
        self.content=self.reader.text
        self.template=temps.si
        self.selects:list[int]=[1,2,3]
        self.sameFile:bool=True
        self.FILENAME:str=''
    
    def write_energy(self):
        engNums=self.reader.read_energys()
        """匹配并保存各种校正能量"""
        for i,eng in enumerate(engNums):
            self.template=self.template.replace(f'<E{i}>',f'{eng:>15.8f}')
            
    def write_coord(self):
        COORDS=[]
        for atom in self.mol.atoms:
            symbol=atom.symbol
            x,y,z=atom.coord
            COORDS.append(f'{symbol:4} {x:14.8f} {y:14.8f} {z:14.8f}')
        self.template=self.template.replace('<COORD>','\n'.join(COORDS))

    def write_freq(self):
        freqObj = re.findall(r'^\s+Frequencies\s--\s+(-?\d+\.\d+.+)$', self.content, flags=re.M)
        if freqObj is not None:
            matchContent=f''
            for freq in freqObj:
                freqs=re.split(r' +',freq)
                freqs=[f'{e:>15}' for e in freqs]
                freqs=''.join(freqs)
                matchContent+=f'{freqs}\n'
            self.template=self.template.replace('<FREQ>',matchContent[:-1])
        else:
            printer.wrong('Match Frequencies Error')

    def build(self):
        """
        导出分子的各种信息到txt文件
        1. 坐标 coord
        2. 能量 energy
        3. 频率 freq
        """
        if 1 in self.selects:
            self.write_coord()
        else:
            self.template=self.template.replace('<COORD>\n','')
        if 2 in self.selects:
            self.write_energy()
        else:
            self.template=self.template.replace('<ENERGY>\n','')
        if 3 in self.selects:
            self.write_freq()
        else:
            self.template=self.template.replace('<FREQ>\n','')
        
    def save(self):
        path_=Path(self.reader.path)
        if self.FILENAME=='':self.FILENAME=path_.stem
        self.template=self.template.replace('<FILENAME>',self.FILENAME)
        if self.sameFile:
            path=f'{path_.parent}/SI.txt'
            mode='a'
            self.template+=('-'*60+'\n')
        else:
            path=f'{path_.parent}/{path_.stem}_SI.txt'
            mode='w'
        with open(path,mode,encoding='utf-8') as f:
            f.write(self.template)

    def onShell(self):
        from pywfn.utils import parse_intList
        print('(1:能量,2:坐标,3:频率):')
        msgs=input('请输入需要保存的信息[*]')
        if msgs!='':self.selects=parse_intList(msgs)
        same=input('是否保存到同一文件内？[y]/n')
        if same=='n':self.sameFile=False
        self.save()