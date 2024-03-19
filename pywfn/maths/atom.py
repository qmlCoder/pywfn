import numpy as np

from pywfn import base

def pola2cart(r,t,p):
    """极坐标转直角坐标，半径，仰角，转角"""
    z=r*np.sin(t)
    x=np.cos(t)*np.cos(p)
    y=np.cos(t)*np.sin(p)
    return x,y,z

def get_direct(atom:"base.Atom"):
    """
    根据最大排斥原则，迭代寻找反应方向
    在原子的外层采点，计算目标值，目标值最小作为初猜
    """
    def get_score(point:np.ndarray): # 衡量当前点好坏的一个函数
        point=point.reshape(1,3)
        angles=-np.dot(point,vectors)+1 # 不是真正的角度，但与真正的角度正相关
        mean=np.mean(angles)
        score=sum([ang**2-(ang-mean)**2 for ang in angles])
        return score

    points=np.zeros(shape=(50,3))
    for t in np.linspace(0,np.pi,5):
        for p in np.linalg(0,2*np.pi,10):
            x,y,z=pola2cart(1,t,p)
            points[t*5+p]=[x,y,z]

    vectors=np.concatenate([around.coord-atom.coord for around in atom.neighbors],axis=0)
    
    scores=np.array([get_score(point) for point in points])
    maxScore=np.max(scores)
    index=np.argmax(scores)
    point0=points[index] # 获取初猜
    arounds=np.array([
        [-1,0],
        [ 1,0],
        [0, 1],
        [0,-1]
    ])
    stepSize=np.pi/10
    for i in range(100):
        points=point0+arounds*stepSize
        scores=np.array([get_score(point) for point in points])
        maxScore_=np.max(scores)
        index=np.argmax(scores)
        point0_=points[index]
        if maxScore_-maxScore>1e-5:
            point0=point0_
            maxScore=maxScore_
        else:
            if stepSize<np.pi/100:
                return point0
            else:
                stepSize/=2