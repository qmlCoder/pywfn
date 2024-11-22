"""
obj文件生成器
"""
import numpy as np

class ObjWriter:
    def __init__(self):
        self.verts:np.ndarray|None=None # 顶点
        self.faces:np.ndarray|None=None # 面
    
    def save(self,path:str):
        assert self.verts is not None,"没有提供顶点"
        assert self.faces is not None,"没有提供面"
        with open(path,'w',encoding='utf-8') as f:
            for i,(x,y,z) in enumerate(self.verts):
                f.write(f'v {x:>10.4f}{y:>10.4f}{z:>10.4f}\n')
            for f1,f2,f3 in self.faces:
                f.write(f'f {f1:>4d} {f2:>4d} {f3:>4d}\n')