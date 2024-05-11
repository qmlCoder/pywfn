"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""
from pywfn.base import Mol
from pywfn.atomprop import atomCharge,atomSpin,atomDirect
from pywfn.atomprop.atomCharge import Chrgs
import numpy as np

class Calculator:
    def __init__(self) -> None:
        self.mols:list[Mol]=[]

    def fukui(self,chrg:Chrgs='mulliken')->np.ndarray:
        """计算所有原子的福井函数[n,2]"""
        assert len(self.mols)==3,"应该有三个分子"
        cals=[atomCharge.Calculator(mol) for mol in self.mols]
        natm=len(self.mols[0].atoms)
        chgs=np.zeros(shape=(natm,3)) # 记录所有原子的电荷
        dchg=np.zeros(shape=(natm,2)) # 记录电荷差值
        for c,cal in enumerate(cals):
            res=cal.calculate(chrg)
            chgs[:,c]=res
        dchg[:,0]=-(chgs[:,0]-chgs[:,1])
        dchg[:,1]=-(chgs[:,1]-chgs[:,2])
        return dchg
    
    def parr(self,chrg:Chrgs='mulliken')->np.ndarray:
        """计算所有原子的parr函数[n,2]"""
        assert len(self.mols)==2,"应该有两个分子"
        cals=[atomSpin.Calculator(mol) for mol in self.mols]
        natm=len(self.mols[0].atoms)
        spins=np.zeros(shape=(natm,3))
        for c,cal in enumerate(cals):
            spins[:,c]=cal.calculate(chrg)
        result=np.zeros(shape=(natm,2))
        result[:,0]=spins[:,0]-spins[:,1]
        result[:,1]=spins[:,1]-spins[:,2]
        return spins
    
    def delEng(self):
        from pywfn.atomprop import atomEnergy
        cals=[atomEnergy.Calculator(mol) for mol in self.mols]
        natm=len(cals[0].mol.atoms)
        engs=np.zeros(shape=(natm,3))
        for c,cal in enumerate(cals):
            res=cal.calculate()
            engs[:,c]=res
        result=np.zeros(shape=(natm,2))
        result[:,0]=engs[:,0]-engs[:,1]
        result[:,1]=engs[:,1]-engs[:,2]
        return engs

    def valence(self)->list[float]:
        """计算原子化合价"""
        assert len(self.mols)==1,"只能算一个分子"
        from pywfn.bondprop import bondOrder
        mol=self.mols[0]
        caler=bondOrder.Calculator(mol)
        orders=caler.mayer()
        orderDict={}
        for a1,a2,order in orders:
            a1,a2=int(a1),int(a2)
            key=f'{a1}-{a2}' if a1<a2 else f'{a2}-{a1}'
            orderDict[key]=order
        result=np.zeros(len(mol.atoms))
        for a,atom in enumerate(mol.atoms):
            a1=atom.idx
            valence=0
            for a2 in atom.neighbors:
                key=f'{a1}-{a2}' if a1<a2 else f'{a2}-{a1}'
                valence+=orderDict[key]
            result[a]=valence
        return result

    def freeValence(self,atms:list[int]):
        """计算指定原子的自由价[d,5](atm,x,y,z,val)"""
        assert len(self.mols)==1,"只能算一个分子"
        from pywfn.bondprop import bondOrder
        mol=self.mols[0]
        caler=bondOrder.Calculator(mol)
        # STAND=1.6494
        STAND=4.0
        result=[]
        for atm1 in atms:
            dirs=None
            orders=[]
            for atm2 in mol.atom(atm1).neighbors:
                bond=[atm1,atm2]
                res=caler.dirMayer([bond]) #[d,4](a1,a2,x,y,z,v)
                if dirs is None:dirs=res[:,2:5]
                vals=res[:,5]
                orders.append(vals)
            orderSum=np.array(orders).sum(axis=0)
            for d in range(len(dirs)):
                x,y,z=dirs[d]
                order=STAND-orderSum[d]
                result.append([atm1,x,y,z,order])
        return np.array(result)

    def dirFukui(self,atms:list[int])->np.ndarray:
        """计算指定原子，指定方向的福井函数[n,6](atm,x,y,z,E,N)"""
        assert len(self.mols)==3,"需要三个分子"
        crgs=[mol.charge for mol in self.mols]
        assert crgs[0]<crgs[1]<crgs[2],"电荷顺序不符"
        cals=[atomCharge.Calculator(mol) for mol in self.mols]
        
        vals:list[np.ndarray]=[]
        dirs=None
        idxs=None
        for cal in cals:
            obts=cal.mol.O_obts
            res=cal.dirCharge('mulliken',obts,atms)
            if dirs is None:dirs=res[:,1:4]
            if idxs is None:idxs=res[:,0]
            vals.append(res[:,4])
        num=len(dirs)
        result=np.zeros(shape=(num,6))
        result[:,0]=idxs
        result[:,1:4]=dirs
        result[:,4]=vals[0]-vals[1]
        result[:,5]=vals[1]-vals[2]
        return result
            
