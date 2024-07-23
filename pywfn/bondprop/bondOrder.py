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

class Calculator:
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    def mayer(self,PM=None,bonds=None)->np.ndarray:
        """计算指定键的mayer键级，可以有多个"""
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
        """带有方向的Mayer键级[d,6](a1,a2,x,y,z,v)"""
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
        """
        计算与指定愿原子相邻的键的束缚键级
        """
        dirCaler=direction.Calculator(self.mol)
        dirs=dirCaler.reaction(atm)
        nebs=self.mol.atom(atm).neighbors
        bonds=[[atm,neb] for neb in  nebs]
        obts=self.mol.O_obts
        result=[]
        
        for d in range(len(dirs)):
            CMp=self.mol.projCM(obts,[atm],[dirs[d]],True,True,akeeps=nebs)
            PMp=CM2PM(CMp,obts,self.mol.oE)
            orders=self.mayer(PM=PMp,bonds=bonds)
            x,y,z=dirs[d]
            for a1,a2,val in orders:
                result.append([a1,a2,x,y,z,val])
        return np.array(result)
    
    def piOrder(self):
        """
        计算pi键级，每一个可能的π键计算出一个π键级
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

    def hmo(self):
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