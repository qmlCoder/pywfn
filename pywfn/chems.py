"""
计算一些简单的物理化学性质，不需要输入文件
"""
import numpy as np
from pywfn.data import consts as cs
from pywfn.utils import printer
# 计算反应活化能 reaction activation energy
def reaActEne(K:float,T:float)->float:
    """
    计算化学反应活化能,KJ/mol
    K:催化常数
    T:反应温度 开尔文
    """
    c = 2 * 1e10 * T
    x = np.log(K / c)
    g = -x * cs.R * T * 1e-3  # kJ/mol
    # g1 = g / 4.184
    return g

# 计算立体选择ee值 Stereo selection ee value
def steSelEE(deR:float,deS:float,T:float):
    """
    计算立体选择性ee值
    """
    RT = cs.R * T
    RS = np.exp((deS - deR) * 4185.8518 / RT)
    ee = (RS - 1) / (RS + 1)
    return ee

# 计算一堆分子能量的玻尔兹曼分布，能量差越大，玻尔兹曼分布越平均
def bezm(engl:np.ndarray):
    """
    计算玻尔兹曼分布
    engs:能量列表
    """
    engs=np.array(engl)
    engs=engs-np.min(engs) # Hartee
    engs=engs*4.5e-18 # J
    k=1.380694e-23 # J/K
    T=298.15 # K
    return np.exp(-engs/(k*T))

def onShell():
    printer.options('化学计算器',{
        '1':'计算反应活化能',
        '2':'计算立体选择ee值',
        '3':'计算玻尔兹曼分布',
    })
    opt=input('请输入对应选项：')
    match opt:
        case '1':
            K=float(input('催化常数K:'))
            T=float(input('反应温度T:'))
            print('反应活化能:',reaActEne(K,T),'kJ/mol')
        case '2':
            deR=float(input('反应热力学参数:'))
            deS=float(input('立体选择热力学参数:'))
            T=float(input('反应温度T:'))
            print('立体选择ee值:',steSelEE(deR,deS,T))
        case '3':
            engs=input('能量列表(空格分隔):').split()
            engs=[float(i) for i in engs]
            engs=np.array(engs)
            print('玻尔兹曼分布:',bezm(engs))
        case _:
            return