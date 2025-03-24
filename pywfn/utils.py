# 将需要的每个模块共同需要的计算都放在这里面
import re
import numpy as np

from pywfn import config


def parse_intList(string,start=0,offset=0)->list[int]:
    """将字符串解析为数字列表"""
    res = []
    for each in re.split(r',|，', string):
        content = each.split('-')
        if len(content) == 1:
            res.append(int(each) - 1 + start + offset)
        else:
            res += [int(i) - 1 + start + offset for i in range(int(content[0]), int(content[1]) + 1)]
    return res

def parse_obtList(string:str,nbas:int):
    """
    将轨道字符转为正常字符
    """
    alps:list[str]=re.findall(r'a\d+',string)
    bets:list[str]=re.findall(r'b\d+',string)
    for alp in alps:
        string=string.replace(alp,alp[1:])
    for bet in bets:
        fbet=f'{int(bet[1:])+nbas}'
        string=string.replace(bet,fbet)
    return parse_intList(string)
    


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
from pywfn import config
import inspect
class Printer:
    def __init__(self) -> None:
        self.ifDebug=config.IF_DEBUG
        self.console=Console()
        self.tables:dict[str,Table]={}

    def __call__(self,text,style='',end='\n'):
        self.console.print(text,style=style,end=end)

    def warn(self,text):
        self.__call__(text,style='#fcc419')

    def info(self,text):
        self.__call__(text,style='#91a7ff')

    def err(self,text):
        self.__call__(text,style='#e03131')
    
    def res(self,text):
        self.__call__(text,style='#a9e34b')
        self.bar()

    def bar(self,len=40):
        self.__call__('-'*len,style='#6741d9')
    
    def log(self,text): # 打印日志
        frame=inspect.stack()[1]
        text=f'{frame.function}|{frame.filename}:{frame.lineno}\n{text}'
        self.__call__(text,style='#228be6')
    
    def print(self,text):
        self.__call__(text)
    
    def multi(self,texts:list[str]):
        for text in texts:
            self.console.print(text,end='')
        self.console.print('')
    
    def vector(self,tip:str,vector:np.ndarray):
        if not config.IF_SHELL:return
        nums=[f'{v:.2f}' for v in vector]
        numStr=','.join(nums)
        self.console.print(tip,end='',style='#dee2e6')
        self.console.print(numStr,style='#f59f00')
    
    def options(self,title,opts:dict[str,str]):
        if not config.IF_SHELL:return
        print(f'{title}')
        print(f"{'选项'}{'功能'}")
        # table=Table(title=title,box=box.SIMPLE_HEAD,title_style="bold black on white")
        # table.add_column('选项',justify="left")
        # table.add_column('功能',justify="left")
        for idx,text in opts.items():
            # table.add_row(f'{idx: >2}. ',f'{text}')
            print(f'{idx: >2}. {text}')
        # print('')
        # self.console.print(table)
    
    def shell(self,text:str): # 以shell方式运行的时候才会打印
        if not config.IF_SHELL:return
        self.console.print(text)

printer=Printer()

class Caler:
    def __init__(self):
        self.showMesg=False

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
    
from typing import Any
def l2i(nums:list[Any]):
    return [int(num) for num in nums]

def l2f(nums:list[Any]):
    return [float(num) for num in nums]

import sys
def err_stop(tip):
    print(tip)
    input("按任意键退出...")
    sys.exit()

def chkArray(array:np.ndarray,shape:list[int|None]):
    """对数组的形状进行检测

    Args:
        array (np.ndarray): 传入数组
        shape (list[int]): 是否符合形状
    """
    if len(shape)!=len(array.shape):return False # 如果维度不同肯定就不符合了
    for i in range(len(shape)):
        if shape[i] is None:continue # 如果是None就不检测
        if shape[i]!=array.shape[i]:return False # 如果形状不同就不符合
    return True