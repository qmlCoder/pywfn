"""
obj文件生成器
"""
import numpy as np

class ObjWriter:
    def __init__(self,verts:np.ndarray,faces:np.ndarray):
        self.verts=verts # 顶点
        self.faces=faces # 面
    
    def save(self,path:str):
        with open(path,'w',encoding='utf-8') as f:
            for i,(x,y,z) in enumerate(self.verts):
                f.write(f'v {x:>10.4f}{y:>10.4f}{z:>10.4f}\n')
            for f1,f2,f3 in self.faces:
                f.write(f'f {f1:>4d} {f2:>4d} {f3:>4d}\n')