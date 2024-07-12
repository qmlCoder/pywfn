"""
显示CDXML文件，将CDXML文件保存为图片
编辑每个键的颜色，原本的颜色就不要使用了，直接上新的颜色
"""

from xml.etree import ElementTree
from xml.etree.ElementTree import Element
from dataclasses import dataclass

@dataclass
class Node:
    nid:int
    x:float
    y:float

@dataclass
class Edge:
    eid:int
    B:int
    E:int
    color:int


class Tool:
    def __init__(self,path:str) -> None:
        self.path=path
        self.tree=ElementTree.parse(path)
        self.root=self.tree.getroot()
        self.colorTable=self.root.find('colortable')
            
    def set_color(self,paras:dict[str:list[float]]):
        """
        设置指定边的颜色，n1,n2,r,g,b
        """
        keys=list(paras.keys())
        edges=self.root.findall('.//b')
        cidx=10
        for edge in edges:
            B=edge.attrib['B']
            E=edge.attrib['E']
            key0=f'{B}-{E}'
            key1=f'{E}-{B}'
            if key0 in keys:
                r,g,b=paras[key0]
                self.colorTable.append(Element('color',{'r':str(r),'g':str(g),'b':str(b)}))
                edge.set('color',str(cidx))
                cidx+=1
            elif key1 in keys:
                r,g,b=paras[key1]
                self.colorTable.append(Element('color',{'r':str(r),'g':str(g),'b':str(b)}))
                edge.set('color',str(cidx))
                cidx+=1
    
    def save(self,name:str):
        self.tree.write(name)