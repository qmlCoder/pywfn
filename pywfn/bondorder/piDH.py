"""
轨道分解法+HMO计算π键级
"""
from pywfn.base import Mol,Atom,Bond
import numpy as np
from pywfn import utils
from pywfn import maths
from pywfn.bondorder import Caler

class Calculator(Caler):
    def __init__(self,mol:Mol) -> None:
        self.mol=mol
        self.bond:list[int]=None

    def get_order(self, center:Atom,around:Atom, orbital:int, direction):  # 计算一个轨道的键级
        '''
        计算两个原子间某个轨道的键级
        中心原子，相邻原子，轨道，方向'''
        if around.symbol=='H':
            return 0
        if np.linalg.norm(direction)==0:
            return 0
        As=self.mol.As[orbital]

        centerPos=center.coord
        aroundPos=around.coord
        bondVector=aroundPos-centerPos
        crossVector=np.cross(utils.normalize(bondVector),utils.normalize(direction))
        

        centerTs=center.OC[:,orbital].copy()
        aroundTs=around.OC[:,orbital].copy()
        if np.linalg.norm(centerTs)==0 or np.linalg.norm(aroundTs)==0:
            return 0
        
        centerPsProj=center.get_pProj(direction, orbital) # 先计算投影
        aroundPsProj=around.get_pProj(direction, orbital)

        nx=bondVector
        ny=crossVector
        nz=direction
        centerRes=utils.coordTrans(nx,ny,nz, centerPsProj)
        aroundRes=utils.coordTrans(nx,ny,nz, aroundPsProj)

        centerPZs=[each[-1] for each in centerRes]
        aroundPZs=[each[-1] for each in aroundRes]

        pOrder=sum([cpz*apz/As for cpz,apz in zip(centerPZs,aroundPZs)])*self.mol.oE
        return pOrder
        # return order

    def get_orders(self,center:Atom,around:Atom,orbitals:list[int],direction):
        """计算每两个原子之间的键级"""
        orders=[self.get_order(center,around,orbital,direction) for orbital in orbitals] # 所有的占据轨道都计算键级
        return orders

    def calculate(self) -> dict:
        """指定一个键，计算该键的键级"""
        idx1,idx2=self.bond
        centerAtom=self.mol.atom(idx1)
        aroundAtom=self.mol.atom(idx2)
        self.centerPIdx=[i for i,l in enumerate(centerAtom.layers) if 'P' in l] #原子的序数
        self.aroundPIdx=[i for i,l in enumerate(aroundAtom.layers) if 'P' in l]
        O_obts=self.mol.O_obts
        normal=centerAtom.get_Normal(aroundAtom) # 原子的法向量
        if len(centerAtom.neighbors)==3: # 如果原子有法向量(sp2)
            orders=self.get_orders(centerAtom,aroundAtom,O_obts,normal)
            return sum(orders)
        else: # 如果没有法向量的话，则需要找轨道方向作为基础方向
            # 从高到低计算轨道方向，遇到有垂直于键轴的则作为轨道方向
            orbitalDirection=None
            for o in O_obts[::-1]:
                orbitalDirection=centerAtom.get_obtWay(o)
                if np.linalg.norm(orbitalDirection)==0:
                    continue
                bondDirection=self.mol.bond(centerAtom.idx, aroundAtom.idx).vector
                if maths.vector_angle(orbitalDirection, bondDirection,trans=True)>0.4: # 夹角要很大
                    break
            if orbitalDirection is not None:
                orders1=self.get_orders(centerAtom,aroundAtom,O_obts,orbitalDirection)
                crossDirection=np.cross(orbitalDirection, bondDirection)
                orders2=self.get_orders(centerAtom,aroundAtom,O_obts,crossDirection)
                return sum(orders1),sum(orders2)
    
    def resStr(self)->str:
        order=self.calculate()
        return f'{order}'
