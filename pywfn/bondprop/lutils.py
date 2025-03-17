"""
公用的工具函数
"""
import numpy as np

from pywfn.base import Atom,Mole
from pywfn.maths import vector_angle
from pywfn.utils import printer
from pywfn.maths.atom import get_sCont
from pywfn.atomprop import direction
from pywfn.gridprop import density

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

def judgeOrbital(mol:Mole,atm1:int,atm2:int,obt:int,dirCaler:direction.Calculator)->int:
    """
    判断一个轨道是否为π轨道，几何方法
    centerAtom,aroundAtom:原子对象
    orbital:分子轨道序数
    normal:原子法向量[其实可以是任意方向]
    """
    
    # print('判断pi轨道:',atm1,atm2)
    atom1=mol.atom(atm1)
    atom2=mol.atom(atm2)
    if atom1.symbol=='H' or atom2.symbol=='H':
        # print(obt+1,'轨道是H原子轨道,无法判断')
        return 0
    # 根据键轴上的电子密度判断
    RA=atom1.coord*0.7+atom2.coord*0.3
    RB=atom1.coord*0.3+atom2.coord*0.7
    densCaler=density.Calculator(mol)
    dens=densCaler.molDens(np.array([RA,RB]),0)
    if np.max(dens)<0.01:return 0
    # 1. 根据s轨道和p轨道的贡献
    sContCenter=get_sCont(mol,atm1,obt)
    sContAround=get_sCont(mol,atm1,obt)
    # print('s轨道贡献:',sContCenter,sContAround)
    if sContCenter>0.01 or sContAround>0.01:
        # print(obt+1,'s轨道贡献:',sContCenter,sContAround)
        return 0
    # 2. p轨道的方向要处在垂直分子平面方向
    cenDir=dirCaler.maxWeave(atm1,obt,'P[XYZ]')
    aroDir=dirCaler.maxWeave(atm2,obt,'P[XYZ]')
    # print('p轨道方向:',atm1,cenDir)
    # print('p轨道方向:',atm2,aroDir)
    
    if cenDir is None and aroDir is None:
        # print(obt+1,'p轨道方向:',cenDir,aroDir)
        return 0
    normal=dirCaler.normal(atm1)
    if cenDir is None: cenDir=normal
    if aroDir is None: aroDir=normal
    centerAngle=vector_angle(cenDir,normal) # type: ignore # 计算分子平面和p轨道方向的夹角
    aroundAngle=vector_angle(aroDir,normal) # type: ignore
    
    if abs(0.5-centerAngle)<0.3 or abs(0.5-aroundAngle)<0.3:
        # print(obt+1,'夹角:',centerAngle,aroundAngle)
        return 0
    # 以上两个条件都满足的可以认为是π轨道
    if (0.5-centerAngle)*(0.5-aroundAngle)>0:
        return 1
    else:
        return -1

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

