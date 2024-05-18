# 将需要的每个模块共同需要的计算都放在这里面
import re
import numpy as np

from pywfn import config


def parse_atmList(string)->list[int]:
    """将字符串解析为数字列表"""
    res = []
    for each in re.split(r',|，', string):
        content = each.split('-')
        if len(content) == 1:
            res.append(int(each) - 1)
        else:
            res += [int(i) - 1 for i in range(int(content[0]), int(content[1]) + 1)]
    return res

def get_vertical(v1,v2):
    '''获得垂直于两向量的单位向量'''
    v=np.cross(v1,v2)
    return v/np.linalg.norm(v)

def get_points_between_two_pos(pos1,pos2,n):
    '''计算两点之间的点的坐标'''
    dp=pos2-pos1
    all_pos=[pos1+dp/(n-1)*i for i in range(n)]
    return np.array(all_pos)


def coordTrans(nx,ny,nz,coords):
    """对一个向量实行空间坐标变换
    nx,ny,nz: 另一组基坐标
    coord:要尽心变换的坐标/向量
    """
    e=np.array([[1,0,0],[0,1,0],[0,0,1]]).T # 空间坐标变换
    nx,ny,nz=normalize(nx),normalize(ny),normalize(nz)
    e_=np.array([nx,ny,nz]).T
    A_=np.dot(np.linalg.inv(e_),e)
    ps__=np.array([np.dot(A_,p) for p in coords])
    return ps__

arounds=np.array([
        [-1,0,0],
        [+1,0,0],
        [0,-1,0],
        [0,+1,0],
        [0,0,-1],
        [0,0,+1]
    ])

def differ_function(posan1,posan2): #计算电子分布差值图
    return (posan1+posan2)**2-(posan1**2+posan2**2)/2


def get_gridPoints(range,step,ball=False): 
    '''获取空间格点[3,n]'''
    points=[]
    for x in np.arange(-range,range,step):
        for y in np.arange(-range,range,step):
            for z in np.arange(-range,range,step):
                if ball:
                    distance=(x**2+y**2+z**2)**0.5
                    if distance>range:
                        continue
                points.append([x,y,z])
    return np.array(points).T

def k_means(des,points):
    '''定义k-means聚类算法函数,返回每一类对应的索引'''
    N=len(des)
    for epoch in range(5):
        types=[[] for n in range(N)]
        idxs=[[] for n in range(N)]
        for _,each in enumerate(points):
            distance=np.linalg.norm((des-each),axis=1) #计算每个点与中心点的距离(范数)，两个数值
            idx=np.argmin(distance) # 获取距离最近的索引
            types[idx].append(each) # 在索引位置添加点（分类）
            idxs[idx].append(_)
        types=[np.array(each) for each in types]
        pos=np.array([np.mean(each,axis=0) if len(each)>0 else des[i] for i,each in enumerate(types)]) # 获取两个类别的中心点坐标
        des=pos  # 更新目标点
    return idxs

def list_remove(l,x):
    new_l=[]
    for each in l:
        if each!=x:
            new_l.append(each)
    return new_l

def connectH(atom,connections):
    '''
    如果连接的有H,就返回C-H键的向量,否则返回None
    '''
    CHVectors=[]
    for each in connections:
        if each['atom_type']=='H':
            CHVectors.append(each['pos']-atom['pos'])
    return CHVectors

def multiple(a,b):
    '''获取两个量之间的倍数关系，必定是大于等于1的'''
    if abs(a)>=abs(b):
        return a/b
    else:
        return b/a

def nodeNum(data):
    '''获取节点曲线中的节点数量，有可能会找多'''
    data[np.where(data==0)]=1e-6
    res=data[:-1]*data[1:]
    return len(np.where(res<=0)[0])

def normalize(vector):
    length=np.linalg.norm(vector)
    if length==0:
        raise
    return vector/length

from rich.console import Console
from rich .table import Table,box
from typing import Sequence,Iterable
from pywfn import config
import inspect
class Printer:
    def __init__(self) -> None:
        self.ifDebug=config.IF_DEBUG
        self.ifShell=config.IF_SHELL
        self.console=Console()
        self.tables:dict[str,Table]={}

    def __call__(self,text,style='',end='\n'):
        if self.ifShell:self.console.print(text,style=style,end=end)

    def warn(self,text):
        self.__call__(text,style='#fcc419')

    def info(self,text):
        self.__call__(text,style='#91a7ff')

    def wrong(self,text):
        self.__call__(text,style='#e03131')
    
    def res(self,text):
        self.__call__(text,style='#a9e34b')
        self.bar()
    
    def print(self,text):
        self.__call__(text)
    
    def multi(self,texts:list[str]):
        self.__call__(text,style='')
        for text in texts:
            self.console.print(text,end='')
        self.console.print('')
    
    def bar(self,len=40):
        self.__call__('-'*len,style='#6741d9')

    def space(self):
        self.__call__('')
    
    def log(self,text):
        if self.ifDebug: # 只有debug模式启动的时候才会打印
            frame=inspect.stack()[1]
            text=f'{frame.function}|{frame.filename}:{frame.lineno}\n{text}'
            self.__call__(text,style='#228be6')
    
    def vector(self,tip:str,vector:np.ndarray):
        if not config.IF_SHELL:return
        vector=[f'{v:.2f}' for v in vector]
        vectorStr=','.join(vector)
        self.console.print(tip,end='',style='#dee2e6')
        self.console.print(vectorStr,style='#f59f00')
    
    def options(self,title,opts:dict[str,str]):
        if not config.IF_SHELL:return
        if title not in self.tables.keys():
            table=Table(title=title,box=box.SIMPLE_HEAD,title_style="bold black on white")
            table.add_column('选项',justify="left")
            table.add_column('功能',justify="left")
            for idx,text in opts.items():
                table.add_row(f'{idx: >2}. ',f'{text}')
            self.tables[title]=table
        table=self.tables[title]
        self.console.print(table)
    
    def track(self,seq:Sequence,tip:str='')->Iterable:
        if config.IF_DEBUG:
            from rich import progress
            return progress.track(seq,description=tip)
        else:
            return seq

printer=Printer()

"""不需要多次实例化的类可以使用类方法"""
import time
class TimeCounter:
    """统计每个函数执行的时间"""
    counts={}
    
    def __init__(self) -> None:
        pass
    
    @classmethod
    def funcTime(cls,func):
        t0=time.time()
        ans=func()
        t1=time.time()
        
        return ans

def vectStr(vector:None|np.ndarray,f=4):
    if vector is None:
        return 'None'
    else:
        resStr=','.join([f'{v:>8.4f}' for v in vector])
        return f'({resStr})'