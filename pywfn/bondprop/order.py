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
from pywfn.bondprop import lutils
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
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        atmuls=self.mol.atoms.atmuls
        bonds=self.mol.bonds.ats
        orders=lutils.mayer(PM,SM,atmuls,bonds)
        # orders=[]
        # for bond in self.mol.bonds:
        #     a1,a2=bond.ats
        #     u1,l1=self.mol.atom(a1).obtBorder
        #     u2,l2=self.mol.atom(a2).obtBorder
        #     vals=OM[u1:l1,u2:l2]
        #     order=np.sum(vals)
        #     orders.append([a1,a2,order])
        # order = np.array(orders)
        return bonds,orders
    
    def multiCenter(self,atms:list[int]):
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
        return result
    
    def dirMayer(self,bond:tuple[int,int],dirs:np.ndarray,aleep=False,lkeep=False)->np.ndarray:
        """计算带有方向的Mayer键级，指定一个键的多个方向

        Args:
            bonds (list[int]): 指定要计算哪个键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        obts=self.mol.O_obts
        result=[]
        a1,a2=bond
        if a1>a2:a1,a2=a2,a1
        for d,dir in enumerate(dirs):
            # CMp=projCM(self.mol,obts,[a1,a2],np.array([dir,dir]),False,True)
            CMp=projCM(self.mol,obts,[a1,a2],np.array([dir,dir]),False,False)
            PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
            orders=lutils.mayer(PMp,self.mol.SM,self.mol.atoms.uls,self.mol.bonds.ats)
            for (a1_,a2_),order in zip(self.mol.bonds.ats,orders):
                if a1!=a1_ or a2!=a2_:continue
                x,y,z=dir
                result.append([a1,a2,x,y,z,order])
        return np.array(result)
    
    def boundMayer(self,atm:int,dirs:np.ndarray)->tuple[list[tuple[int,int]],np.ndarray]:
        """计算与指定原子相邻的键的束缚键级，束缚键级投影指定原子，保留邻接原子，清除其它原子

        Args:
            atm (int): 指定原子编号

        Returns:
            np.ndarray: 束缚键级`[a1,a2,x,y,z,val]`
        """
        if config.SHOW_LEVEL>=1:
            printer.info(f'计算原子{atm}的束缚键级')
        nebs=self.mol.atom(atm).neighbors
        obts=self.mol.O_obts
        result=[]
        bonds=[]
        for d in range(len(dirs)):
            CMp=projCM(self.mol,obts,[atm],dirs[d,np.newaxis],True,True,akeeps=nebs)
            PMp=CM2PM(CMp,obts,self.mol.oE)
            orders=lutils.mayer(PMp,self.mol.SM,self.mol.atoms.uls,self.mol.bonds.ats)
            x,y,z=dirs[d]
            for (a1,a2),val in zip(self.mol.bonds.ats,orders):
                if int(a1) not in nebs+[atm]:continue
                if int(a2) not in nebs+[atm]:continue
                bonds.append((a1,a2))
                result.append(val)
                if config.SHOW_LEVEL>=1:
                    print(f"{a1}{a2}{val}")
        result=np.array(result)
        # print('束缚键级',atm,result)
        return bonds,result
    
    def pi_pocv(self)->tuple[list[tuple[int,int]],np.ndarray]:
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

    def pi_smo(self,bond:list[int])->np.ndarray:
        """根据分子轨道挑选法计算pi键级

        Args:
            atm1 (int): 键的第一个原子
            atm2 (int): 键的第二个原子

        Returns:
            np.ndarray: pi键级
        """
        dirCaler=direction.Calculator(self.mol)
        atm1,atm2=bond
        atom1=self.mol.atom(atm1)
        atom2=self.mol.atom(atm2)
        
        # normal=dirCaler.normal(atm1)
        # assert normal is not None,f'原子{atm1}没有法向量'
        obts=self.mol.O_obts
        CMs=np.zeros_like(self.mol.CM) # 拷贝一份，然后将不是π轨道的那些变成0

        piObts=[]
        for atom in [atom1,atom2]: # 修改每个原子对应的系数矩阵
            if atom.symbol=='H':continue
            u,l=atom.obtBorder
            for obt in obts:
                judgeRes=lutils.judgeOrbital(self.mol,atm1,atm2,obt,dirCaler)
                if judgeRes==0:continue # 如果不是π轨道，将判断的不是pi轨道的轨道系数置为0
                CMs[u:l,obt]=self.mol.CM[u:l,obt]
                if obt not in piObts:piObts.append(obt)
        print(f'挑选的pi轨道有：{piObts}')
        oe=self.mol.oE
        PMs=lutils.CM2PM(CMs,obts,oe)
        SM=self.mol.SM
        PS=PMs@SM

        u1,l1=atom1.obtBorder
        u2,l2=atom2.obtBorder
        order=np.sum(PS[u1:l1,u2:l2]*PS[u2:l2,u1:l1])
        return order

    def hmo(self)->np.ndarray:
        """计算休克尔分子轨道法的键级

        Returns:
            np.ndarray: HMO键级
        """
        atms=self.mol.heavyAtoms
        natm=len(atms)
        BM,es,CM,occs=hmo(self.mol)
        # 3.构建键级矩阵
        nmat=CM.shape[0]
        CM[:,nmat//2:]=0 
        orders=[]
        for i,ai in enumerate(atms):
            for j,aj in enumerate(atms):
                if i>=j:continue
                if BM[i,j]==0:continue
                dist=self.mol.DM[ai-1,aj-1]
                if dist>1.7*1.8897:continue
                order=np.sum(CM[i,:]*CM[j,:])*2
                # print(i,j,CM[i,:],CM[j,:],order)
                orders.append([ai,aj,order])
        orders=np.array(orders)
        return np.abs(orders)

    
    
    # 分解键级
    def decompose(self,bond:list[int],dobt:int=-1):
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