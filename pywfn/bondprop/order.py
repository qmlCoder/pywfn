"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mole
from pywfn.atomprop import direction
from pywfn.maths import CM2PM
from pywfn.maths.mol import hmo,projCM
from pywfn.utils import printer
from pywfn.cli import Shell
from pywfn.bondprop import lutils

import numpy as np
from itertools import product
from pywfn import config
from collections import defaultdict
from itertools import product

class Calculator:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol

    # mayer键级
    def mayer(self)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算mayer键级，mayer键级是基础键级，很多方法的键级都是基于mayer键级计算出来的

        Args:
            PM (np.ndarray, optional): 密度矩阵，如不指定则使用分子默认密度矩阵
            bonds (list[list[int]], optional): 可指定要计算的键，若不指定则使用所有键

        Returns:
            np.ndarray: 所有键的键级，形状为:[d,3](a1,a2,order)
        """
        # 获取密度矩阵 P
        PM=self.mol.PM
        # 获取重叠矩阵 S
        SM=self.mol.SM
        atmuls=self.mol.atoms.atmuls
        bonds=self.mol.bonds.ats
        orders=lutils.mayer(PM,SM,atmuls,bonds)
        return bonds,orders
    
    def MCBO(self,atms:list[int])->float:
        """计算多中心键级"""
        PS=self.mol.PM@self.mol.SM
        uls=[range(*self.mol.atom(atm).obtBorder)for atm in atms]
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
    
    def dirOrder(self,bond:tuple[int,int],dir_:np.ndarray,aleep=False,lkeep=False)->float:
        """计算带有方向的Mayer键级，指定一个键的多个方向

        Args:
            bonds (list[int]): 指定要计算哪个键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        obts=self.mol.O_obts
        a1,a2=bond
        CMp=projCM(self.mol,obts,[a1,a2],np.array([dir_,dir_]),False,False)
        PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
        orders=lutils.mayer(PMp,self.mol.SM,self.mol.atoms.uls,self.mol.bonds.ats)
        for (a1_,a2_),order in zip(self.mol.bonds.ats,orders):
            if a1!=a1_ or a2!=a2_:continue
            return order.item()
        raise ValueError(f"没有找到键级，bond={bond},dir={dir_}")
    
    def boundOrder(self,atm:int,dir_:np.ndarray)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算与指定原子相邻的键的束缚键级，束缚键级投影指定原子，保留邻接原子，清除其它原子

        Args:
            atm (int): 指定原子编号
            dir_ (np.ndarray): 指定原子投影的方向，形状为[3,](x,y,z)

        Returns:
            np.ndarray: 束缚键级`[a1,a2,x,y,z,val]`
        """
        if config.SHOW_LEVEL>=1:
            printer.info(f'计算原子{atm}的束缚键级')
        nebs=self.mol.atom(atm).neighbors
        obts=self.mol.O_obts
        result=[]
        bonds=[]
        CMp=projCM(self.mol,obts,[atm],dir_.reshape(1,3),True,True,akeeps=nebs)
        PMp=CM2PM(CMp,obts,self.mol.oE)
        orders=lutils.mayer(PMp,self.mol.SM,self.mol.atoms.uls,self.mol.bonds.ats)
        for (a1,a2),val in zip(self.mol.bonds.ats,orders):
            if int(a1) not in nebs+[atm]:continue
            if int(a2) not in nebs+[atm]:continue
            bonds.append((a1,a2))
            result.append(val)
            if config.SHOW_LEVEL>=1:
                print(f"{a1}{a2}{val}")
        result=np.array(result)
        return bonds,result
    
    def piOrder_pocv(self)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算分子的所有pi键键级，每个可能的pi键计算出一个键级

        Returns:
            np.ndarray: 返回ndarray数组
        """
        dirCaler=direction.Calculator(self.mol)
        atms=[]
        dirs=[]
        for atom in self.mol.atoms:
            normal=dirCaler.normal(atom.idx) # 原子的法向量
            if normal is None:continue
            atms.append(atom.idx)
            dirs.append(normal)
        dirs=np.array(dirs)
        PMp=projCM(self.mol,self.mol.O_obts,atms,dirs,False,False)
        PMp=CM2PM(PMp,self.mol.O_obts,self.mol.oE)
        
        bonds=self.mol.bonds.ats
        orders=lutils.mayer(PMp,self.mol.SM,self.mol.atoms.uls,self.mol.bonds.ats)
        orders[orders<0]=0
        orders=np.sqrt(orders)
        return bonds,orders

    def piOrder_smo(self,bond:tuple[int,int])->float:
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
        dirCaler=direction.Calculator(self.mol)
        densCaler=density.Calculator(self.mol)
        atm1,atm2=bond
        atom1=self.mol.atom(atm1)
        atom2=self.mol.atom(atm2)
        if atom1.symbol=='H' or atom2.symbol=='H':
            return 0
        # 根据键轴上的电子密度判断
        RA=atom1.coord*0.7+atom2.coord*0.3
        RB=atom1.coord*0.3+atom2.coord*0.7
        dens=densCaler.molDens(np.array([RA,RB]),0)[0]
        if np.max(dens)<0.01:
            return 0
        obts=self.mol.O_obts
        CMs=np.zeros_like(self.mol.CM) # 拷贝一份，然后将不是π轨道的那些变成0

        piObts=[]
        for atom in [atom1,atom2]: # 修改每个原子对应的系数矩阵
            u,l=atom.obtBorder
        # u1,l1=atom1.obtBorder
        # u2,l2=atom1.obtBorder
            for obt in obts:
                # judgeRes=lutils.judgeOrbital(self.mol,atm1,atm2,obt,dirCaler,densCaler)
                
                # 1. 根据s轨道和p轨道的贡献
                sContCenter=get_sCont(self.mol,atm1,obt)
                sContAround=get_sCont(self.mol,atm1,obt)
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
                normal=dirCaler.normal(atm1)
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
                CMs[u:l,obt]=self.mol.CM[u:l,obt].copy()
                piObts.append(obt)
        print(f'挑选的pi轨道有：{set(piObts)}')
        oe=self.mol.oE
        PMs=lutils.CM2PM(CMs,obts,oe)
        SM=self.mol.SM
        result=lutils.mayer(PMs,SM,self.mol.atoms.atmuls,[bond])
        return result[0].item()

    def HMO(self)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算休克尔分子轨道法的键级

        Returns:
            np.ndarray: HMO键级
        """
        atms=self.mol.heavyAtoms
        natm=len(atms)
        BM,es,CM,occs=hmo(self.mol)
        nmat=CM.shape[0]
        CM[:,nmat//2:]=0 
        orders=[]
        bonds=[]
        for i,ai in enumerate(atms):
            for j,aj in enumerate(atms):
                if i>=j:continue
                if BM[i,j]==0:continue
                dist=self.mol.DM[ai-1,aj-1]
                if dist>1.7*1.8897:continue
                order=np.sum(CM[i,:]*CM[j,:])*2
                bonds.append((ai,aj))
                orders.append(order)
        orders=np.array(orders)
        return bonds,np.abs(orders)

    
    
    # 分解键级
    def decompose(self,bond:tuple[int,int],dobt:int=-1):
        """
        键级分解，将两个原子的轨道分解到指定的局部坐标系中，然后根据每种键的重叠模式计算键级
        将原子轨道基函数的系数按照角动量进行分组
        atm1,atm2:组成键的两个原子
        keeps:每个角动量保留第几个系数，例如{1:[2]}代表p轨道只保留pz
        """
        from pywfn.orbtprop.decom import decomOrbitals
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol) # 方向计算器
        # 计算出两个坐标系
        atm1,atm2=bond
        T1=dirCaler.coordSystem(atm1,atm2) # 两个原子的局部坐标系
        T2=dirCaler.coordSystem(atm2,atm1)
        Ts=[T1,T2]
        
        nmat=self.mol.CM.shape[0]
        if dobt==-1:
            keepList=[
                [[1],[1,0,0],[1,1,1,0,0,0]],
                [[0],[0,0,1],[0,0,0,0,1,0]],
                [[0],[0,1,0],[0,0,0,1,0,0]],
                [[0],[0,0,0],[0,0,0,0,0,1]],
            ]
            
        elif dobt==0:
            keepList=[
                [[1],[1,0,0],[0,0,0,0,0,0]],
                [[0],[0,0,1],[0,0,0,0,0,0]],
                [[0],[0,1,0],[0,0,0,0,0,0]],
                [[0],[0,0,0],[0,0,0,0,0,0]],
            ]
        elif dobt==1:
            keepList=[
                [[1],[1,0,0],[1,1,1,1,1,1]],
                [[0],[0,0,1],[1,1,1,1,1,1]],
                [[0],[0,1,0],[1,1,1,1,1,1]],
                [[0],[0,0,0],[1,1,1,1,1,1]],
            ]
        else:
            raise ValueError("dobt must be 0 or 1 or -1")

        orders=[]
        CMs=[]
        for k,keeps in enumerate(keepList):
            CMt=np.zeros_like(self.mol.CM) # 变换矩阵初始化
            for o in self.mol.O_obts:
                coefDict=defaultdict(list) # 系数字典
                for i in range(nmat):
                    iatm=self.mol.atoAtms[i]
                    ishl=self.mol.atoShls[i]
                    iang=self.mol.atoAngs[i]
                    key=(iatm,ishl,iang)
                    coefDict[key].append(self.mol.CM[i,o])

                for key,val in coefDict.items():
                    iatm,ishl,iang=key
                    rcoefs=np.array(val) # 原始系数
                    if iatm in bond:
                        T=Ts[bond.index(iatm)]
                        tcoefs=decomOrbitals(T,rcoefs,keeps[iang],dtype='bond')
                        # print(rcoefs,tcoefs)
                    else:
                        tcoefs=rcoefs
                    assert len(rcoefs)==len(tcoefs),"长度对不上"
                    
                    coefDict[key]=tcoefs.tolist() # type: ignore
                values=list(coefDict.values())
                CMt[:,o]=np.concatenate(values)
            CMs.append(CMt)
            
            PMt=CM2PM(CMt,self.mol.O_obts,self.mol.oE) # 变换的密度矩阵
            bonds=self.mol.bonds.ats
            orders=lutils.mayer(PMt,self.mol.PM,self.mol.atoms.uls,bonds)
            orders[orders<0]=0
            orders=orders**0.5
            result=[]
            for (a1,a2),val in zip(bonds,orders):
                if a1 not in bond:continue
                if a2 not in bond:continue
                result.append(val)
                break
        return np.array(result)