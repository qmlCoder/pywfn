"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol
from pywfn.atomprop import direction
from pywfn.maths import CM2PM
from pywfn.utils import printer
import numpy as np
from itertools import product
from pywfn import config
from collections import defaultdict

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def mayer(self,PM:np.ndarray|None=None,bonds:list[list[int]]|None=None)->np.ndarray:
        """计算mayer键级

        Args:
            PM (np.ndarray, optional): 密度矩阵，如不指定则使用分子默认密度矩阵
            bonds (list[list[int]], optional): 可指定要计算的键，若不指定则使用所有键

        Returns:
            np.ndarray: 指定键的键级
        """
        # 获取密度矩阵 P
        if PM is None:PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        orders=[]
        if bonds is None:
            bonds=[bond.ats for bond in self.mol.bonds]
        for a1,a2 in bonds:
            u1,l1=self.mol.atom(a1).obtBorder
            u2,l2=self.mol.atom(a2).obtBorder
            vals=OM[u1:l1,u2:l2]
            order=np.sum(vals)
            orders.append([a1,a2,order])
        order = np.array(orders)
        # print(order)
        return order
    
    def dirMayer(self,bonds:list[list[int]])->np.ndarray:
        """计算带有方向的Mayer键级

        Args:
            bonds (list[list[int]]): 通过二维整数列表指定要计算哪些键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        dirCaler=direction.Calculator(self.mol)
        obts=self.mol.O_obts
        result=[]
        for a1,a2 in bonds:
            dirs=dirCaler.reaction(a1)
            atms=[a1]*len(dirs)
            if a1>a2:a1,a2=a2,a1
            for d in range(len(dirs)):
                CMp=self.mol.projCM(obts,[a1,a2],[dirs[d],dirs[d]],False,True)
                PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
                a1_,a2_,order=self.mayer(PM=PMp,bonds=[[a1,a2]]).flatten()
                assert a1==a1_ and a2==a2_,"原子不对应"
                x,y,z=dirs[d]
                result.append([a1,a2,x,y,z,order])
        return np.array(result)
    
    def boundMayer(self,atm:int)->np.ndarray:
        """计算与指定原子相邻的键的束缚键级

        Args:
            atm (int): 指定原子编号

        Returns:
            np.ndarray: 束缚键级
        """
        dirCaler=direction.Calculator(self.mol)
        dirs=dirCaler.reaction(atm)
        nebs=self.mol.atom(atm).neighbors
        bonds=[[atm,neb] for neb in  nebs]
        obts=self.mol.O_obts
        result=[]
        
        for d in range(len(dirs)):
            CMp=self.mol.projCM(obts,[atm],[dirs[d]],False,True,akeeps=nebs)
            PMp=CM2PM(CMp,obts,self.mol.oE)
            orders=self.mayer(PM=PMp,bonds=bonds)
            x,y,z=dirs[d]
            for a1,a2,val in orders:
                result.append([a1,a2,x,y,z,val])
        return np.array(result)
    
    def piOrder(self)->np.ndarray:
        """计算pi键键级，每个可能的pi键计算出一个键级

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
        
        PMp=self.mol.projCM(self.mol.O_obts,atms,dirs,False,False)
        PMp=CM2PM(PMp,self.mol.O_obts,self.mol.oE)
        result=self.mayer(PM=PMp)
        orders=result[:,-1]
        orders[orders<0]=0
        orders=np.sqrt(orders)
        result[:,-1]=orders
        return result

    def hmo(self)->np.ndarray:
        """计算休克尔分子轨道法的键级

        Returns:
            np.ndarray: HMO键级
        """
        # 1.建立键连矩阵
        atms=self.mol.heavyAtoms
        natm=len(atms)
        BM=np.zeros(shape=(natm,natm)) # 键连矩阵
        DM=np.zeros_like(BM) # 键长矩阵
        for i,j in product(range(natm),range(natm)):
            a1,a2=atms[i],atms[j]
            if a1>=a2:continue
            bond=self.mol.atom(a1).coord-self.mol.atom(a2).coord
            dist=np.linalg.norm(bond)
            DM[i,j]=dist
            DM[j,i]=dist
            if dist>1.7*1.889:continue
            BM[i,j]=1.0
            BM[j,i]=1.0
        # 2.求解
        e,C=np.linalg.eig(BM) # 矩阵对角化
        nele=int(len(atms)-self.mol.charge) #电子数量
        idxs=np.argsort(e)[:nele//2] # 占据轨道
        CM=C[:,idxs] # 每一列对应一个特征向量
        
        # 3.构建键级矩阵
        result=[]
        OM=np.zeros_like(BM)
        for i,j in product(range(natm),range(natm)):
            order=np.sum(CM[i,:]*CM[j,:])*2
            OM[i,j]=order
            if i>=j:continue
            if DM[i,j]>1.7*1.889:continue
            result.append([atms[i],atms[j],order])
        result=np.array(result)
        return np.abs(result)

    def multiCenter(self,atms:list[int]):
        """计算多中心键级"""
        pass
    
    # 分解键级
    def decomOrder(self,atms:list[int],keeps:dict[int,list[int]]):
        """
        键级分解，将两个原子的轨道分解到指定的局部坐标系中，然后根据每种键的重叠模式计算键级
        将原子轨道基函数的系数按照角动量进行分组
        atm1,atm2:组成键的两个原子
        keeps:每个角动量保留第几个系数，例如{1:[2]}代表p轨道只保留pz
        """
        import matplotlib.pyplot as plt
        from pywfn.atomprop import direction
        dirCaler=direction.Calculator(self.mol) # 方向计算器
        # 计算出两个坐标系
        atm1,atm2=atms
        T1=dirCaler.coordSystem(atm1,atm2) # 两个原子的局部坐标系
        T2=dirCaler.coordSystem(atm2,atm1)
        Ts=(T1,T2)
        # CM=self.mol.CM.copy()
        CM=np.zeros_like(self.mol.CM)
        nmat=CM.shape[0]

        for o in self.mol.O_obts:
            coefDict=defaultdict(list) #系数字典
            for i in range(nmat):
                iatm=self.mol.obtAtms[i]
                ishl=self.mol.obtShls[i]
                iang=self.mol.obtAngs[i]
                key=(iatm,ishl,iang)
                if iatm in atms:
                    coefDict[key].append(self.mol.CM[i,o])
                else:
                    coefDict[key].append(0)

            for key,val in coefDict.items():
                iatm,ishl,iang=key
                rcoefs=np.array(val)
                if iatm in atms:
                    tcoefs=decomOrbitals(Ts[atms.index(iatm)],rcoefs,keeps[iang])
                    # print(iatm,ishl,iang,rcoefs,tcoefs)
                else:
                    tcoefs=rcoefs
                assert len(rcoefs)==len(tcoefs),"长度对不上"
                
                coefDict[key]=tcoefs
            values=list(coefDict.values())
            CM[:,o]=np.concatenate(values)
        # print(CM)
        # fig,axs=plt.subplots(1,2)
        # axs[0].matshow(self.mol.CM,vmin=-1,vmax=1)
        # axs[1].matshow(CM,vmin=-1,vmax=1)
        # plt.show()
        PM=CM2PM(CM,self.mol.O_obts,self.mol.oE)
        orders=self.mayer(PM)
        # orders[:,-1]=np.sqrt(orders[:,-1])
        return orders

    def onShell(self):
        while True:
            printer.options('键级计算',{
                '1':'mayer键级',
                '2':'方向mayer键级',
                '3':'pi键级',
                '4':'HMO键级',
            })
            opt=input('请输入序号选择要计算的键级：')
            if opt=='1':
                orders=self.mayer()
                for a1,a2,val in orders:
                    print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
            elif opt=='2':
                opt=input('请输入需要计算的键，例如(1-2): ')
                a1,a2=opt.split('-')
                result=self.dirMayer(bonds=[[int(a1),int(a2)]])
                for a1,a2,x,y,z,val in result:
                    print(f'{int(a1):>2d}-{int(a2):>2d}({x:>8.4f} {y:>8.4f} {z:>8.4f}):{val:>8.4f}')
            elif opt=='3':
                orders=self.piOrder()
                for a1,a2,val in orders:
                    print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
            elif opt=='4':
                orders=self.hmo()
                for a1,a2,val in orders:
                    print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
            elif opt=='':
                break
            else:
                printer.warn('无效选项!')


def gtf(cords:np.ndarray,lmn:tuple[int,int,int])->np.ndarray:
    """计算高斯型波函数

    Args:
        cords (np.ndarray): 格点坐标
        lmn (tuple[int,int,int]): 角动量分量

    Returns:
        np.ndarray: 波函数数值
    """
    x=cords[0,:]
    y=cords[1,:]
    z=cords[2,:]
    l,m,n=lmn
    r2=x**2+y**2+z**2
    alp=2.0
    facs=[1,1,3]
    fac=facs[l]*facs[m]*facs[n]
    ang=l+m+n
    Nm=(2*alp/np.pi)**(3/4)*np.sqrt((4*alp)**ang/fac)
    val=x**l * y**m * z**n * np.exp(-alp*r2)*Nm
    return val


def decomOrbitals(T:np.ndarray,coefs:np.ndarray,keeps:list):
    if len(coefs)==1:
        return decomOrbitalS(T,coefs,keeps)
    elif len(coefs)==3:
        return decomOrbitalP(T,coefs,keeps)
    elif len(coefs)==6:
        return decomOrbitalD(T,coefs,keeps)
    else:
        printer.warn("不支持的数组长度")
        return coefs

def decomOrbitalS(T:np.ndarray,coefs:np.ndarray,keeps:list):
    if keeps:
        return coefs
    else:
        return np.array([0.])

# 分解P轨道
def decomOrbitalP(T:np.ndarray,rcoefs:np.ndarray,keeps:list):
    npos=3
    cords=np.random.rand(3,npos) #随机生成6个点
    zs1=np.zeros(shape=(npos,3))
    zs1[:,0]=gtf(cords,(1,0,0))
    zs1[:,1]=gtf(cords,(0,1,0))
    zs1[:,2]=gtf(cords,(0,0,1))

    zs2=np.zeros(shape=(npos,3))
    zs2[:,0]=gtf(T@cords,(1,0,0))
    zs2[:,1]=gtf(T@cords,(0,1,0))
    zs2[:,2]=gtf(T@cords,(0,0,1))

    Mr=np.linalg.inv(zs2)@zs1
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    for i in range(3):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs
    return fcoefs

# 分解D轨道
def decomOrbitalD(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    npos=6
    cords=np.random.rand(3,npos) #随机生成6个点
    zs1=np.zeros(shape=(npos,6))
    zs1[:,0]=gtf(cords,(2,0,0))
    zs1[:,1]=gtf(cords,(0,2,0))
    zs1[:,2]=gtf(cords,(0,0,2))
    zs1[:,3]=gtf(cords,(1,1,0))
    zs1[:,4]=gtf(cords,(1,0,1))
    zs1[:,5]=gtf(cords,(0,1,1))

    zs2=np.zeros(shape=(npos,6))
    zs2[:,0]=gtf(T@cords,(2,0,0))
    zs2[:,1]=gtf(T@cords,(0,2,0))
    zs2[:,2]=gtf(T@cords,(0,0,2))
    zs2[:,3]=gtf(T@cords,(1,1,0))
    zs2[:,4]=gtf(T@cords,(1,0,1))
    zs2[:,5]=gtf(T@cords,(0,1,1))

    Mr=np.linalg.inv(zs2)@zs1
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    for i in range(6):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs
    return fcoefs