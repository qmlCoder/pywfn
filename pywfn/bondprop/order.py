"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base.mole import Mole
from pywfn.atomprop import direction
from pywfn.maths import CM2PM
from pywfn.maths.mol import hmo,projCM
from pywfn.utils import printer
from pywfn.bondprop import lutils

import numpy as np
from itertools import product
from pywfn import config
from collections import defaultdict
from itertools import product
from pywfn.utils import chkArray
from pywfn import core
from pywfn.moleprop.orbital import Deco

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.caler=core.bondprop.order.Calculator(mole.mole) # type: ignore # 核心计算器

    # mayer键级
    def mayer(self)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算mayer键级，mayer键级是基础键级，很多方法的键级都是基于mayer键级计算出来的

        Args:
            PM (np.ndarray, optional): 密度矩阵，如不指定则使用分子默认密度矩阵
            bonds (list[list[int]], optional): 可指定要计算的键，若不指定则使用所有键

        Returns:
            np.ndarray: 所有键的键级，形状为:[d,3](a1,a2,order)
        """
        return self.caler.mayer()
    
    def MCBO(self,atms:list[int])->float:
        """计算多中心键级"""
        PS=self.mole.PM@self.mole.SM
        uls=[range(*self.mole.atom(atm).obtBorder)for atm in atms]
        result=0.0
        for each in product(*uls):
            val=1.0
            for i in range(len(each)):
                if i!=len(each)-1:
                    val*=PS[each[i],each[i+1]]
                else:
                    val*=PS[each[i],each[0]]
            result+=val
        unit=abs(result)/result
        result=unit*abs(result)**(1/len(atms))
        return result.item()
    
    def pocv(self,dirs:dict[int,list[float]],keep_other_atm:bool,keep_other_sym:bool)->float:
        """计算带有方向的Mayer键级，指定一个键的多个方向

        Args:
            bonds (list[int]): 指定要计算哪个键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        return self.caler.pocv(dirs,keep_other_atm,keep_other_sym,"mayer") # type: ignore
    
    def deco(self,decos:dict[int,Deco],ctype:str='mulliken'):
        """分解键级"""
        return self.caler.deco(decos,ctype)
    
    def bound(self,atm:int,dir:list[float])->tuple[list[tuple[int,int]],np.ndarray]:
        """计算与指定原子相邻的键的束缚键级，束缚键级投影指定原子，保留邻接原子，清除其它原子

        Args:
            atm (int): 指定原子编号
            dir (list[float]): 指定原子投影的方向，形状为[3,](x,y,z)

        Returns:
            np.ndarray: 束缚键级`[a1,a2,x,y,z,val]`
        """
        return self.caler.bound(atm,dir)
    
    def pi_pocv(self) -> tuple[dict[int,list[float]], np.ndarray]:
        """计算分子的所有pi键键级，每个可能的pi键计算出一个键级
        """
        dirs,vals=self.caler.pi_pocv()
        return  dirs,np.array(vals)
    
    def pi_deco(self):
        deco,vals=self.caler.pi_deco()
        return  deco,np.array(vals) # type: ignore

    def pi_smo(self,bond:tuple[int,int])->float:
        """根据分子轨道挑选法计算pi键级

        Args:
            atm1 (int): 键的第一个原子
            atm2 (int): 键的第二个原子

        Returns:
            np.ndarray: pi键级
        """
        from pywfn.maths.atom import get_sCont
        from pywfn.maths import vector_angle
        from pywfn.gridprop import density
        dirCaler=direction.Calculator(self.mole)
        densCaler=density.Calculator(self.mole)
        atm1,atm2=bond
        atom1=self.mole.atom(atm1)
        atom2=self.mole.atom(atm2)
        if atom1.symbol=='H' or atom2.symbol=='H':
            return 0
        # 根据键轴上的电子密度判断
        RA=atom1.coord*0.7+atom2.coord*0.3
        RB=atom1.coord*0.3+atom2.coord*0.7
        dens=densCaler.molDens(np.array([RA,RB]),0)[0]
        if np.max(dens)<0.01:
            return 0
        obts=self.mole.O_obts
        CMs=np.zeros_like(self.mole.CM) # 拷贝一份，然后将不是π轨道的那些变成0

        piObts=[]
        for atom in [atom1,atom2]: # 修改每个原子对应的系数矩阵
            u,l=atom.obtBorder
        # u1,l1=atom1.obtBorder
        # u2,l2=atom1.obtBorder
            for obt in obts:
                # judgeRes=lutils.judgeOrbital(self.mol,atm1,atm2,obt,dirCaler,densCaler)
                
                # 1. 根据s轨道和p轨道的贡献
                sContCenter=get_sCont(self.mole,atm1,obt)
                sContAround=get_sCont(self.mole,atm1,obt)
                # print('s轨道贡献:',sContCenter,sContAround)
                if sContCenter>0.01 or sContAround>0.01:
                    continue
                # 2. p轨道的方向要处在垂直分子平面方向
                cenDir=dirCaler.maxWeave(atm1,obt,'P[XYZ]')
                aroDir=dirCaler.maxWeave(atm2,obt,'P[XYZ]')
                # print('p轨道方向:',atm1,cenDir)
                # print('p轨道方向:',atm2,aroDir)
                
                if cenDir is None and aroDir is None:
                    continue
                normal=dirCaler.normal_vector(atm1)
                if cenDir is None: cenDir=normal
                if aroDir is None: aroDir=normal
                centerAngle=vector_angle(cenDir,normal) # type: ignore # 计算分子平面和p轨道方向的夹角
                aroundAngle=vector_angle(aroDir,normal) # type: ignore
                
                if abs(0.5-centerAngle)<0.3 or abs(0.5-aroundAngle)<0.3:
                    continue
                # 以上两个条件都满足的可以认为是π轨道
                # if (0.5-centerAngle)*(0.5-aroundAngle)>0:
                #     return 1
                # else:
                #     return -1
                # if judgeRes==0:continue # 如果不是π轨道，将判断的不是pi轨道的轨道系数置为0
                

                # CMs[u1:l1,obt]=self.mol.CM[u1:l1,obt]
                # CMs[u2:l2,obt]=self.mol.CM[u2:l2,obt]
                CMs[u:l,obt]=self.mole.CM[u:l,obt].copy()
                piObts.append(obt)
        print(f'挑选的pi轨道有：{set(piObts)}')
        oe=self.mole.oE
        PMs=lutils.CM2PM(CMs,obts,oe)
        SM=self.mole.SM
        result=lutils.mayer(PMs,SM,self.mole.atoms.atmuls,[bond])
        return result[0].item()

    def HMO(self)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算休克尔分子轨道法的键级

        Returns:
            np.ndarray: HMO键级
        """
        atms=self.mole.heavyAtoms
        natm=len(atms)
        BM,es,CM,occs=hmo(self.mole)
        nmat=CM.shape[0]
        CM[:,nmat//2:]=0 
        orders=[]
        bonds=[]
        for i,ai in enumerate(atms):
            for j,aj in enumerate(atms):
                if i>=j:continue
                if BM[i,j]==0:continue
                dist=self.mole.DM[ai-1,aj-1]
                if dist>1.7*1.8897:continue
                order=np.sum(CM[i,:]*CM[j,:])*2
                bonds.append((ai,aj))
                orders.append(order)
        orders=np.array(orders)
        return bonds,np.abs(orders)