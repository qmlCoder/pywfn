"""
显示CDXML文件，将CDXML文件保存为图片
需要用到matplotlib了呀，本来不想用的呢
"""

import matplotlib.pyplot as plt
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
        self.colors=list(self.root.find('colortable')) # type: ignore # 所有的颜色
        self.page=self.root.find('page')
        self.frags=self.page.findall('fragment') # type: ignore
        self.nodes:list[Node]=[]
        self.edges:list[Edge]=[]
        self.find_data()
    
    def find_data(self):
        for frag in self.frags:
            self.parse_frag(frag)

    def parse_frag(self,frag:Element):
        for child in frag:
            if child.tag=='n':
                nid=child.attrib['id']
                x,y=child.attrib['p'].split(' ')
                self.nodes.append(Node(int(nid),float(x),-float(y)))
            if child.tag=='b':
                eid=child.attrib['id']
                B=child.attrib['B']
                E=child.attrib['E']
                if 'color' in child.attrib.keys():
                    color=int(child.attrib['color'])-2
                else:
                    color=1
                self.edges.append(Edge(int(eid),int(B),int(E),color))

    def get_node(self,nid:int):
        for node in self.nodes:
            if node.nid==nid:
                return node
        
    def get_color(self,idx:int):
        attrib=self.colors[idx].attrib
        r,g,b=[float(attrib[c]) for c in ('r','g','b')]
        return r,g,b

            
    def draw(self):
        # 先画边
        for edge in self.edges:
            B=edge.B
            E=edge.E
            NB=self.get_node(B)
            NE=self.get_node(E)
            color=self.get_color(edge.color)
            assert NB is not None,"没有找到开始节点"
            assert NE is not None,"没有找到结束节点"
            plt.plot([NB.x,NE.x],[NB.y,NE.y],color=color)
        
        # 再画节点
        xs=[]
        ys=[]
        for node in self.nodes:
            xs.append(node.x)
            ys.append(node.y)
            plt.text(node.x,node.y,str(node.nid),color='red')
        plt.scatter(xs,ys,s=10)
        plt.show()