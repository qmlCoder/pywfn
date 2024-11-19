"""
键级也不止一种，都在这里实现吧
"""
from pywfn.base import Mol
from pywfn.atomprop import direction
from pywfn.maths import CM2PM
from pywfn.maths.mol import hmo,projCM
from pywfn.utils import printer

import numpy as np
from itertools import product
from pywfn import config
from collections import defaultdict
from pywfn.bondprop import lutils

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    # mayer键级
    def mayer(self,PM:np.ndarray|None=None)->np.ndarray:
        """计算mayer键级，mayer键级是基础键级，很多方法的键级都是基于mayer键级计算出来的

        Args:
            PM (np.ndarray, optional): 密度矩阵，如不指定则使用分子默认密度矩阵
            bonds (list[list[int]], optional): 可指定要计算的键，若不指定则使用所有键

        Returns:
            np.ndarray: 所有键的键级，形状为:[d,3](a1,a2,order)
        """
        # 获取密度矩阵 P
        if PM is None:PM=self.mol.PM
        # 获取重叠矩阵
        SM=self.mol.SM
        PS=PM@SM
        OM=PS*PS.T
        orders=[]
        for bond in self.mol.bonds:
            a1,a2=bond.ats
            u1,l1=self.mol.atom(a1).obtBorder
            u2,l2=self.mol.atom(a2).obtBorder
            vals=OM[u1:l1,u2:l2]
            order=np.sum(vals)
            orders.append([a1,a2,order])
        order = np.array(orders)
        # print(order)
        return order
    
    def dirMayer(self,bond:list[int])->np.ndarray:
        """计算带有方向的Mayer键级

        Args:
            bonds (list[int]): 指定要计算哪个键级

        Returns:
            np.ndarray: 返回数组形状为:[d,6](a1,a2,x,y,z,v),其中d为键级方向数量
        """
        dirCaler=direction.Calculator(self.mol)
        obts=self.mol.O_obts
        result=[]
        a1,a2=bond
        dirs=dirCaler.reactions(a1) #a1的反应方向
        if a1>a2:a1,a2=a2,a1
        for d,dir in enumerate(dirs):
            CMp=projCM(self.mol,obts,[a1,a2],[dir,dir],False,True)
            PMp=CM2PM(CMp,self.mol.O_obts,self.mol.oE)
            for a1_,a2_,order in self.mayer(PM=PMp):
                if a1!=a1_ or a2!=a2_:continue
                x,y,z=dir
                result.append([a1,a2,x,y,z,order])
        return np.array(result)
    
    def boundMayer(self,atm:int,dirs:np.ndarray)->np.ndarray:
        """计算与指定原子相邻的键的束缚键级，束缚键级投影指定原子，保留邻接原子，清除其它原子

        Args:
            atm (int): 指定原子编号

        Returns:
            np.ndarray: 束缚键级`[a1,a2,x,y,z,val]`
        """
        # if dirs is None:
        #     dirCaler=direction.Calculator(self.mol)
        #     dirs=dirCaler.reaction(atm) # 获取原子的反应方向
        nebs=self.mol.atom(atm).neighbors
        # bonds=[[atm,neb] for neb in  nebs]
        obts=self.mol.O_obts
        result=[]
        
        for d in range(len(dirs)):
            # CMp=projCM(self.mol,obts,[atm],[dirs[d]],False,True,akeeps=nebs)
            CMp=projCM(self.mol,obts,[atm],[dirs[d]],True,True,akeeps=nebs)
            PMp=CM2PM(CMp,obts,self.mol.oE)
            orders=self.mayer(PM=PMp)
            x,y,z=dirs[d]
            for a1,a2,val in orders:
                if int(a1) not in nebs+[atm]:continue
                if int(a2) not in nebs+[atm]:continue
                # print(a1,a2,x,y,z,val)
                result.append([a1,a2,x,y,z,val])
        result=np.array(result)
        return result
    
    def pi_pocv(self)->np.ndarray:
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
        
        PMp=projCM(self.mol,self.mol.O_obts,atms,dirs,False,False)
        PMp=CM2PM(PMp,self.mol.O_obts,self.mol.oE)
        result=self.mayer(PM=PMp)
        orders=result[:,-1]
        orders[orders<0]=0
        orders=np.sqrt(orders)
        result[:,-1]=orders
        return result

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
        
        normal=dirCaler.normal(atm1)

        obts=self.mol.O_obts
        CMs=np.zeros_like(self.mol.CM) # 拷贝一份，然后将不是π轨道的那些变成0

        piObts=[]
        for atom in [atom1,atom2]: # 修改每个原子对应的系数矩阵
            if atom.symbol=='H':continue
            a_1,a_2=atom.obtBorder
            for obt in obts:
                judgeRes=lutils.judgeOrbital(self.mol,atm1,atm2,obt,normal)
                if judgeRes==0:continue # 如果是π轨道
                CMs[a_1:a_2,obt]=self.mol.CM[a_1:a_2,obt]
                if obt not in piObts:piObts.append(obt)
        print(f'挑选的pi轨道有：{piObts}')
        oe=self.mol.oE
        PMs=lutils.CM2PM(CMs,obts,oe)
        SM=self.mol.SM
        PS=PMs@SM
        
        a1,a2=atom1.obtBorder
        b1,b2=atom2.obtBorder
        order=np.sum(PS[a1:a2,b1:b2]*PS[b1:b2,a1:a2])
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

    def multiCenter(self,bond:list[int]):
        """计算多中心键级"""
        pass
    
    # 分解键级
    def decompose(self,bond:tuple[int,int]):
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
        atm1,atm2=bond
        T1=dirCaler.coordSystem(atm1,atm2) # 两个原子的局部坐标系
        T2=dirCaler.coordSystem(atm2,atm1)
        T2=T1
        Ts=(T1.T,T2.T)

        CMt=np.zeros_like(self.mol.CM) # 变换矩阵初始化
        nmat=CMt.shape[0]
        sig_keeps={
            0:[0], #s
            1:[1], #px,py,pz
            2:[]  #xx,yy,zz,xy,xz,yz
        } # 每个角动量保持的索引

        sig_keeps={
            0:[0], #s
            1:[1], #px,py,pz
            2:[1]  #xx,yy,zz,xy,xz,yz
        } # 每个角动量保持的索引

        piz_keeps={
            0:[], 
            1:[2], 
            2:[2,5]  #xx,yy,zz,xy,xz,yz
        } # 每个角动量保持的索引

        pix_keeps={
            0:[], 
            1:[0], 
            2:[0,3]  #xx,yy,zz,xy,xz,yz
        } # 每个角动量保持的索引
        det_keeps={
            0:[], #s
            1:[], #px,py,pz
            2:[4] #xx,yy,zz,xy,xz,yz
        } # 每个角动量保持的索引
        keepList=[sig_keeps,piz_keeps,pix_keeps,det_keeps]
        nameList=['sigma','pi_z','pi_x','delta']
        # keeps=pi_skeeps
        orders=[]
        for k,keeps in enumerate(keepList):
            # if k!=1:continue
            for o in self.mol.O_obts:
                coefDict=defaultdict(list) #系数字典
                for i in range(nmat):
                    iatm=self.mol.obtAtms[i]
                    ishl=self.mol.obtShls[i]
                    iang=self.mol.obtAngs[i]
                    key=(iatm,ishl,iang)
                    if iatm in bond:
                        coefDict[key].append(self.mol.CM[i,o])
                    else:
                        coefDict[key].append(0)

                for key,val in coefDict.items():
                    iatm,ishl,iang=key
                    rcoefs=np.array(val)
                    if iatm in bond:
                        tcoefs=decomOrbitals(Ts[bond.index(iatm)],rcoefs,keeps[iang])
                    else:
                        tcoefs=rcoefs
                    assert len(rcoefs)==len(tcoefs),"长度对不上"
                    
                    coefDict[key]=tcoefs
                values=list(coefDict.values())
                CMt[:,o]=np.concatenate(values)
            PMt=CM2PM(CMt,self.mol.O_obts,self.mol.oE) # 变换的密度矩阵
            results=self.mayer(PMt)
            for a1,a2,val in results:
                if a1 not in bond:continue
                if a2 not in bond:continue
                order=val.item()**0.5
            orders.append(order)
        return orders

    def onShell(self):
        while True:
            printer.options('键级计算',{
                '1':'mayer键级',
                '2':'方向mayer键级',
                '3':'pi键级(POCV)',
                '4':'pi键级(SMO)',
                '5':'HMO键级',
                '6':'分解键级'
            })
            opt=input('选择要计算的键级：')
            match opt:
                case '1':
                    orders=self.mayer()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '2':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        a1,a2=opt.split('-')
                        result=self.dirMayer(bond=[int(a1),int(a2)])
                        for a1,a2,x,y,z,val in result:
                            print(f'{int(a1):>2d}-{int(a2):>2d}({x:>8.4f} {y:>8.4f} {z:>8.4f}):{val:>8.4f}')
                case '3':
                    orders=self.pi_pocv()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '4':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        a1,a2=opt.split('-')
                        order=self.pi_smo(bond=[int(a1),int(a2)])
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{order:>8.4f}')
                case '5':
                    orders=self.hmo()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '6':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        a1,a2=opt.split('-')
                        orders=self.decompose(bond=[int(a1),int(a2)])
                        sig,piz,pix,det=orders
                        print(f'σ : {sig:>10.4f}')
                        print(f'πz: {piz:>10.4f}')
                        print(f'πx: {pix:>10.4f}')
                        print(f'δ : {det:>10.4f}')
                case _:
                    break


def gtf(cords,l,m,n)->np.ndarray:
    x=cords[0,:]
    y=cords[1,:]
    z=cords[2,:]
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
        return coefs

def decomOrbitalS(T:np.ndarray,coefs:np.ndarray,keeps:list[int]):
    if keeps:
        return coefs
    else:
        return np.array([0.])

# 分解P轨道
def decomOrbitalP(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int])->np.ndarray:
    """分解P轨道

    Args:
        T (np.ndarray): 基坐标，每一行代表一个方向
        rcoefs (np.ndarray): 原始函数空间的基函数系数
        keeps (list[int]): 保留的角动量

    Returns:
        np.ndarray: 分解之后的轨道系数
    """
    npos=3
    cords=np.random.rand(3,npos) #随机生成3个点
    zs1=np.zeros(shape=(npos,3))
    zs1[:,0]=gtf(cords,1,0,0)
    zs1[:,1]=gtf(cords,0,1,0)
    zs1[:,2]=gtf(cords,0,0,1)

    zs2=np.zeros(shape=(npos,3))
    zs2[:,0]=gtf(T@cords,1,0,0)
    zs2[:,1]=gtf(T@cords,0,1,0)
    zs2[:,2]=gtf(T@cords,0,0,1)

    # Mt=(np.linalg.inv(zs2)@zs1)

    Mr=np.linalg.inv(zs2)@zs1
    Mi=np.linalg.inv(Mr)
    # print(rcoefs)
    tcoefs=Mr@rcoefs # 根据函数空间基组1下的系数获取函数空间基组2下的系数
    # print(tcoefs)
    for i in range(3):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs # 根据修改后的函数空间基组2下的系数得到函数空间基组1下的系数
    # print(tcoefs)
    # print(fcoefs)
    return fcoefs

# 分解D轨道
def decomOrbitalD(T:np.ndarray,rcoefs:np.ndarray,keeps:list[int]):
    npos=6
    cords=np.random.rand(3,npos) #随机生成6个点
    zs1=np.zeros(shape=(npos,6))
    zs1[:,0]=gtf(cords,2,0,0)
    zs1[:,1]=gtf(cords,0,2,0)
    zs1[:,2]=gtf(cords,0,0,2)
    zs1[:,3]=gtf(cords,1,1,0)
    zs1[:,4]=gtf(cords,1,0,1)
    zs1[:,5]=gtf(cords,0,1,1)

    zs2=np.zeros(shape=(npos,6))
    zs2[:,0]=gtf(T@cords,2,0,0)
    zs2[:,1]=gtf(T@cords,0,2,0)
    zs2[:,2]=gtf(T@cords,0,0,2)
    zs2[:,3]=gtf(T@cords,1,1,0)
    zs2[:,4]=gtf(T@cords,1,0,1)
    zs2[:,5]=gtf(T@cords,0,1,1)

    Mr=np.linalg.inv(zs2)@zs1
    Mi=np.linalg.inv(Mr)
    tcoefs=Mr@rcoefs
    for i in range(6):
        if i in keeps:continue
        tcoefs[i]=0.0
    fcoefs=Mi@tcoefs
    return fcoefs