"""
公用的工具函数
"""
import numpy as np

from pywfn.base import Atom,Mol
from pywfn.maths import vector_angle
from pywfn.utils import printer
from pywfn.maths.atom import get_sCont
from pywfn.atomprop import direction

def CM2PM(CM,orbital:list[int],oe:int)->np.ndarray:
    """
    根据系数矩阵构建密度矩阵
    CM:系数矩阵,如果是开壳层的话,列数是行数的两倍[n,n]/[n,2n]
    n:分子轨道占据电子数
    """
    PMs=CM2PMs(CM,orbital,oe)
    return np.sum(PMs,axis=0)

def CM2PMs(CM,orbital:list[int],oe:int):
    """
    构建三维密度矩阵，不要空轨道，形状为[占据轨道数,原子轨道数,原子轨道数]
    """
    A=(CM[:,orbital].T)[:,:,np.newaxis]
    B=(CM[:,orbital].T)[:,np.newaxis,:]
    return A@B*oe #用矩阵乘法的形式直接构建矩阵可比逐元素计算快多了

def judgeOrbital(mol:Mol,atm1:int,atm2:int,obt:int,normal:np.ndarray)->int:
    """
    判断一个轨道是否为π轨道,几何方法
    centerAtom,aroundAtom:原子对象
    orbital:分子轨道序数
    normal:原子法向量[其实可以是任意方向]
    """
    atom1=mol.atom(atm1)
    atom2=mol.atom(atm2)
    if atom1.symbol=='H' or atom2.symbol=='H':
        return 0
    # 1. 根据s轨道和p轨道的贡献
    sContCenter=get_sCont(mol,atm1,obt)
    sContAround=get_sCont(mol,atm1,obt)
    if sContCenter>0.01 or sContAround>0.01:
        return 0
    # 2. p轨道的方向要处在垂直分子平面方向
    dirCaler=direction.Calculator(mol)
    cenDir=dirCaler.maxWeave(atm1,obt,'p')
    aroDir=dirCaler.maxWeave(atm2,obt,'p')

    if cenDir is None or aroDir is None:
        return 0
    centerAngle=vector_angle(cenDir,normal)
    aroundAngle=vector_angle(aroDir,normal)
    if abs(0.5-centerAngle)<0.3 or abs(0.5-aroundAngle)<0.3:
        return 0
    # 以上两个条件都满足的可以认为是π轨道
    if (0.5-centerAngle)*(0.5-aroundAngle)>0:
        return 1
    else:
        return -1

def judgeOrbital_(centerAtom:Atom,aroundAtom:Atom,orbital:int,normal)-> int:
    """判断某一个键的分子轨道是不是π轨道,成键返回1，反键返回-1，否则返回0"""
    

    centerTs=centerAtom.pLayersTs(orbital)
    aroundTs=aroundAtom.pLayersTs(orbital)
    centerPs=[np.array(centerTs[i:i+3]) for i in range(0,len(centerTs),3)]
    aroundPs=[np.array(aroundTs[i:i+3]) for i in range(0,len(aroundTs),3)]
    normal=centerAtom.get_Normal(aroundAtom)
    centerPs_=centerAtom.get_pProj(normal, orbital) # 先计算投影
    aroundPs_=aroundAtom.get_pProj(normal, orbital)
    centerPs,centerPs_=np.array(centerPs),np.array(centerPs_) # 将系数转为数组
    aroundPs,aroundPs_=np.array(aroundPs),np.array(aroundPs_)
    centerPs,centerPs_=centerPs.sum(axis=0),centerPs_.sum(axis=0) # 将不同组的p求和
    aroundPs,aroundPs_=aroundPs.sum(axis=0),aroundPs_.sum(axis=0)
    centerL,centerL_=np.linalg.norm(centerPs),np.linalg.norm(centerPs_) #p轨道系数组成向量的长度
    aroundL,aroundL_=np.linalg.norm(aroundPs),np.linalg.norm(aroundPs_)
    # centerRatio=np.linalg.norm(centerPs_)/np.linalg.norm(centerPs) # 投影后与投影前的比例
    # aroundRatio=np.linalg.norm(aroundPs_)/np.linalg.norm(aroundPs) 
    centerRatio=np.divide(centerL_,centerL,out=np.zeros_like(centerL),where=centerL!=0)
    aroundRatio=np.divide(aroundL_,aroundL,out=np.zeros_like(aroundL),where=aroundL!=0)

    
    centerScont=centerAtom.get_sCont(orbital)
    aroundScont=aroundAtom.get_sCont(orbital)
    if centerScont>0.001 or aroundScont>0.001:
        return 0 # s贡献太大的不是

    if centerRatio<=0.1 or aroundRatio<=0.1:
        return 0
    if vector_angle(centerPs_,aroundPs_)<=0.5:
        return 1
    else:
        return -1

def printOrders(orders,orbitals):
    """将轨道根据大小排序后输出"""
    orders_=[f'{order:.4f}' for order in orders]
    sortedRes=sorted(list(zip(orbitals,orders_)),key=lambda e:abs(float(e[1])),reverse=True)
    sortedRes=[e for e in sortedRes if abs(float(e[1]))>=0.01]
    sortedRes=list(zip(*sortedRes))    
    formPrint(sortedRes,8,10)

def formPrint(contents:list[list[str]],eachLength:int,lineNum:int=10):
    """
    格式化打印列表内容，contents是一个列表，其中的每一项是一个包含字符串的列表，每个字符串列表长度必须相同
    contents: 待打印的内容，可以想象成一个表
    eachLength: 打印一行中每一项的字符长度
    lineNum: 一行中有多少项
    """
    logs:list[list[str]]=[]
    for content in contents:
        logs.append([])
        for i in range(0,len(content),lineNum):
            text=''.join([f'{each}'.rjust(eachLength,' ') for each in content[i:i+lineNum]])
            logs[-1].append(text)
    if len(logs)==0:
        return
    for i in range(len(logs[0])):
        for log in logs:
            printer.warn(log[i])

