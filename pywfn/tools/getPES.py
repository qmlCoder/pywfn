"""
本脚本用以生成势能面
"""
from pathlib import Path

from pywfn.base import Mol
from pywfn.reader import get_reader,LogReader
from pywfn.data import temps
from collections import namedtuple


Node=namedtuple('Node','id x y')
Bond=namedtuple('Bond','id B E Display')
Text=namedtuple('Text','id t x y')


class Tool:
    def __init__(self,paths:list[Path]) -> None:
        self.paths=paths
        self.route:list[int]
        self.temp=temps.pes
        self.nodes:list[Node]=[]
        self.bonds:list[Bond]=[]
        self.texts:list[Text]=[]
        self.NBS=''
        self.TXS=''
        self.idx=4

    def create(self):
        engs=[]
        for i,idx in enumerate(self.route):
            path=self.paths[idx]
            names,engs=LogReader(str(path)).read_energys()
            eng=engs[-1]
            print(f'{path.name}:{eng=}')
            engs.append(eng*627.51)
        eng0=engs[0]
        engs=[eng-eng0 for eng in engs]
        print(f'转换后{engs=}')
        for i,eng in enumerate(engs):
            self.add_node((i*2)*30+50,300-eng*10)
            self.add_node((i*2+1)*30+50,300-eng*10)

        for i in range(0,len(self.nodes),2):
            nb=self.nodes[i]
            ne=self.nodes[i+1]
            x=(nb.x+ne.x)/2-5
            y=nb.y-5
            self.add_bond(nb.id,ne.id,display='Bold')
            eng=engs[i//2]
            self.add_text(eng,x,y)

        for i in range(1,len(self.nodes)-1,2):
            b=self.nodes[i].id
            e=self.nodes[i+1].id
            self.add_bond(b,e,display='Dash')
        self.weite()

            
    def add_node(self,x,y):
        
        self.nodes.append(Node(id=self.idx,x=x,y=y))
        self.idx+=1
        return self.idx

    def add_bond(self,b:int,e:int,display:str):
        """
        添加键
        """
        self.bonds.append(Bond(id=self.idx,B=b,E=e,Display=display))
        self.idx+=1
        return self.idx
    
    def add_text(self,text,x,y):
        self.texts.append(Text(id=self.idx,t=text,x=x,y=y))
        self.idx+=1
        
    
    def weite(self):
        for node in self.nodes:
            self.NBS+=f'<n id="{node.id}" p="{node.x} {node.y}" Z="1" AS="N"/>\n'
        for bond in self.bonds:
            self.NBS+=f'<b id="{bond.id}" B="{bond.B}" E="{bond.E}" Z="2" Display="{bond.Display}" BS="N"/>\n'
        for text in self.texts:
            self.TXS+=f'<t id="{text.id}" p="{text.x} {text.y}" Z="3" LineHeight="auto"><s font="5" size="10" color="0">{text.t:.1f}</s></t>\n'
        self.temp=self.temp.replace('[NBS]',self.NBS)
        self.temp=self.temp.replace('[TXS]',self.TXS)
        self.temp=self.temp.replace('[BOUND]','0 0 10 20')
        (Path.cwd()/'pes.cdxml').write_text(self.temp)
        # print(self.temp)

    