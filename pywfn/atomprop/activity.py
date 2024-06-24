"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""
from pywfn.base import Mol
from pywfn.atomprop import charge, direction, spin
from pywfn.atomprop.charge import Chrgs
import numpy as np
from pywfn.shell import Shell
from pywfn.maths import CM2PM

class Calculator:
    def __init__(self) -> None:
        self.mols:list[Mol]=[]

    def fukui(self,chrg:Chrgs='mulliken')->np.ndarray:
        """计算所有原子的福井函数[n,2]"""
        assert len(self.mols)==3,"应该有三个分子"
        cals=[charge.Calculator(mol) for mol in self.mols]
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
        cals=[spin.Calculator(mol) for mol in self.mols]
        natm=len(self.mols[0].atoms)
        spins=np.zeros(shape=(natm,3))
        for c,cal in enumerate(cals):
            spins[:,c]=cal.calculate(chrg)
        result=np.zeros(shape=(natm,2))
        result[:,0]=spins[:,0]-spins[:,1]
        result[:,1]=spins[:,1]-spins[:,2]
        return spins
    
    def engDiff(self):
        from pywfn.atomprop import energy
        cals=[energy.Calculator(mol) for mol in self.mols]
        natm=len(cals[0].mol.atoms)
        engs=np.zeros(shape=(natm,3))
        for c,cal in enumerate(cals):
            res=cal.atmEngs()
            engs[:,c]=res
        result=np.zeros(shape=(natm,2))
        result[:,0]=engs[:,0]-engs[:,1]
        result[:,1]=engs[:,1]-engs[:,2]
        return engs

    def valence(self)->np.ndarray:
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

    def freeValence(self,atm:int):
        """计算指定原子的自由价[d,5](atm,x,y,z,val)"""
        assert len(self.mols)==1,"只能算一个分子"
        from pywfn.bondprop import bondOrder
        from pywfn.atomprop import direction
        mol=self.mols[0]
        caler=bondOrder.Calculator(mol)
        # STAND=1.6494
        STAND=3.0
        result=[]
        bmays=caler.boundMayer(atm)
        nebs=mol.atom(atm).neighbors
        for i in range(0,len(bmays),len(nebs)):
            orders=bmays[i:i+len(nebs),-1]
            x,y,z=bmays[i,2:5]
            valence=STAND-sum(orders)
            result.append([atm,x,y,z,valence])
        return np.array(result)

    def delFreeValence(self,atm:int):
        """计算自由价之差"""
        pass

    def dirFukui(self,atms:list[int])->np.ndarray:
        """计算指定原子，指定方向的福井函数[n,6](atm,x,y,z,E,N)"""
        assert len(self.mols)==3,"需要三个分子"
        crgs=[mol.charge for mol in self.mols]
        assert crgs[0]<crgs[1]<crgs[2],"电荷顺序不符"
        cals=[charge.Calculator(mol) for mol in self.mols]
        
        vals:list[np.ndarray]=[]
        dirs=None
        idxs=None
        for cal in cals:
            res=cal.dirCharge('mulliken',atms)
            if dirs is None:dirs=res[:,1:4]
            if idxs is None:idxs=res[:,0]
            vals.append(res[:,4])
        assert dirs is not None,"方向计算出错"
        num=len(dirs)
        result=np.zeros(shape=(num,6))
        result[:,0]=idxs
        result[:,1:4]=dirs
        result[:,4]=vals[0]-vals[1]
        result[:,5]=vals[1]-vals[2]
        return result
    
    def onShell(self,shell:Shell):
        from pywfn.utils import printer
        while True:
            printer.options('原子活性',{
                '1':'福井函数',
                '2':'parr函数',
                '3':'原子能差',
                '4':'化合价',
                '5':'方向自由价',
                '6':'方向福井函数',
            })
            opt=input('选择计算活性类型:')
            chrgMap={'':'mulliken','1':'mulliken','2':'lowdin'}
            if opt=='1':
                self.mols=shell.input.Moles(num=3)
                opt=input('请输入电荷类型：1.mulliken; 2.lowdin')
                if opt not in chrgMap.keys():continue
                result=self.fukui()
                for i,(e,n) in enumerate(result):
                    print(f'{i+1:>3d} {e:>8.4f}{n:>8.4f}')
            elif opt=='2':
                self.mols=shell.input.Moles(num=3)
                opt=input('请输入电荷类型：1.mulliken; 2.lowdin')
                if opt not in chrgMap.keys():continue
                result=self.parr()
                for i,(e,n) in enumerate(result):
                    print(f'{i+1:>3d} {e:>8.4f}{n:>8.4f}')
            elif opt=='3':
                self.mols=shell.input.Moles(num=3)
                result=self.engDiff()
                for i,val in enumerate(result):
                    print(f'{i+1:>3d} {val:>8.4f}')
            elif opt=='4': # 化合价
                self.mols=shell.input.Moles(num=1)
                result=self.valence()
                for i,val in enumerate(result):
                    print(f'{i+1:>3d} {val:>8.4f}')
            elif opt=='5': # 自由价
                self.mols=shell.input.Moles(num=1)
                atms=shell.input.Number(dtype=int,tip='输入原子编号:')
                assert atms is not None,"输入错误"
                for atm in atms:
                    result=self.freeValence(int(atm))
                    for a,x,y,z,v in result:
                        print(f'{a:>3d} {x:>8.4f} {y:>8.4f} {z:>8.4f} {v:>8.4f}')
            elif opt=='6': # 方向fukui函数
                self.mols=shell.input.Moles(num=3)
                atms=shell.input.Number(tip='输入原子编号: ',dtype=int)
                assert atms is not None,"输入错误"
                atms=[int(atm) for atm in atms]
                result=self.dirFukui(atms)
                for i,val in enumerate(result):
                    print(f'{i+1:>3d} {val:>8.4f}')
            else:
                break