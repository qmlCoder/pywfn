from pywfn.base import Mole
from pywfn.reader import LogReader,AnyReader
from pywfn.data import temps
from pywfn.data.elements import elements
from pywfn.base import Geome

from pathlib import Path
import re
import numpy as np


class Tool:
    def __init__(self) -> None:
        pass

    def get_si(self,paths:list[str]) -> str:
        text=''
        for path in paths:
            mol=Mole(LogReader(path))
            assert isinstance(mol.reader,LogReader),"必须是log文件"
            reader:LogReader=mol.reader
            temp=temps.si
            # 生成坐标
            COORDS=[]
            syms=mol.geome.syms
            xyzs=mol.geome.xyzs
            xyzs=xyzs.copy()*0.529177
            for sym,(x,y,z) in zip(syms,xyzs):
                COORDS.append(f'{sym:<12}{x:16.8f}{y:16.8f}{z:16.8f}')
            temp=temp.replace('<COORD>','\n'.join(COORDS))

            # 生成能量
            engNums=reader.read_energys()
            """匹配并保存各种校正能量"""
            for i,eng in enumerate(engNums):
                temp=temp.replace(f'<E{i}>',f'{eng:>15.8f}')

            # 生成频率
            freqs=reader.read_freqs()
            freqStr=''
            for i,freq in enumerate(freqs):
                freqStr+=f'{freq:>20.2f}'
                if (i+1)%3==0 and (i+1)!=len(freqs):freqStr+='\n'
            temp=temp.replace('<FREQ>',freqStr)
            # 生成名称
            temp=temp.replace('<NAME>',reader.fname)
            text+=f"{temp}\n"
        return text
    
    def split_scan(self,path:str):
        start=0
        read=False
        geomes:list[Geome]=[]
        engs:list[float]=[]
        with open(path,'r') as file:
            for l,line in enumerate(file):
                line=line.rstrip()
                # print(line)
                if line.strip()=='Input orientation:':
                    read=True
                    start=l+5
                    geomes.append(Geome())
                if l>=start and line.strip()=='-'*69:
                    read=False
                if read and l>=start:
                    match=re.match(r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)",line)
                    assert match is not None,"匹配结果不正确"
                    atomic,x,y,z=match.groups()
                    sym=elements[int(atomic)].sym
                    x=float(x)
                    y=float(y)
                    z=float(z)
                    # print(sym,x,y,z)
                    geomes[-1].addAtom(sym,np.array([x,y,z]))
                if line[:23]==' SCF Done:  E(RB3LYP) =':
                    # print(line)
                    eng=re.search(r'-?\d+.\d+',line)
                    assert eng is not None,"未匹配到能量"
                    eng=eng.group()
                    engs.append(float(eng))
                if line==' Optimization completed.':
                    break
        moles:list[Mole]=[]
        for g,geome in enumerate(geomes):
            reader=AnyReader()
            reader.geome=geome
            mol=Mole(reader)
            mol._props['energy']=engs[g]
            moles.append(mol)
        return moles

    def split_irc(self):
        pass

    def split_link(self):
        pass
            