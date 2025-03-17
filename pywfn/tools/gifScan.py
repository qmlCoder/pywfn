"""
生成刚性扫描的结构
"""
from pywfn.base import Mole

from pywfn.editor import Editor
from pywfn.data import temps

class Tool:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol
        self.editor=Editor(mol)
        self.temp=temps.gjf
        self.chk="chk"
        self.title="b3lyp/6-31g"
        self.charge=0
        self.spin=1
        self.coord=""
        self.gjfs=[]

    def scan(self,atm1:int,atm2:int,step:int,size:float):
        
        for i in range(step):
            coordl=[]
            coords=self.editor.rotate_bond(atm1,atm2,i*size)/1.889
            for i,(x,y,z) in enumerate(coords):
                s=self.mol.atoms.symbols[i]
                coordl.append(f'{s:>2}{x:>14.8f}{y:>14.8f}{z:>14.8f}')
            gjfText=self.build('\n'.join(coordl))
            self.gjfs.append(gjfText)
            self.temp=temps.gjf

    def build(self,coordStr):
        self.temp=self.temp.replace('<CHK>',self.chk)
        self.temp=self.temp.replace('<TITLE>',self.title)
        self.temp=self.temp.replace('<CHARGE>',f'{self.charge}')
        self.temp=self.temp.replace('<SPIN>',f'{self.spin}')
        self.temp=self.temp.replace('<COORD>',coordStr)
        return self.temp
    
    def save(self,path:str):
        with open(path,'w') as f:
            f.write('--Link1--\n\n'.join(self.gjfs))
        