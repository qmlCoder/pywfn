"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""
from pywfn.base import Mole
from pywfn.atomprop import charge, direction, spin
from pywfn.atomprop.charge import Chrgs
import numpy as np
from pywfn.cli import Shell
from pywfn.maths import CM2PM
from pywfn import config

class Calculator:
    def __init__(self,mol:Mole) -> None:
        self.mol=mol

    # 福井函数
    def fukui(self,molN:Mole,molP:Mole,ctype:str)->np.ndarray:
        """计算所有原子的福井函数

        Args:
            molN (Mole): 负电荷分子
            molP (Mole): 正电荷分子
            ctype (Chrgs, optional): 电荷类型. Defaults to 'mulliken'.

        Returns:
            np.ndarray: 福井函数[n,7](N,N+1,N-1,f-,f+,f0,df)
        """
        mols=[self.mol,molN,molP]
        cals=[charge.Calculator(mol) for mol in mols]
        natm=self.mol.atoms.natm
        vals=np.zeros(shape=(natm,7)) # 记录所有原子的电荷
        for c,cal in enumerate(cals):
            cal.form='charge'
            match ctype:
                case 'mulliken':
                    vals[:,c]=cal.mulliken()
                case 'lowdin':
                    vals[:,c]=cal.lowdin()
                case 'hirshfeld':
                    vals[:,c]=cal.hirshfeld()
                case _:
                    raise ValueError(f'未知的电荷类型{ctype}')
        q0=vals[:,0]
        qn=vals[:,1]
        qp=vals[:,2]
        fn=qp-q0
        fp=q0-qn
        f0=(qp-qn)/2
        df=fp-fn
        vals[:,3]=fn
        vals[:,4]=fp
        vals[:,5]=f0
        vals[:,6]=df
        return vals
    
    # parr函数
    def parr(self,molN:Mole,molP:Mole,ctype:str)->np.ndarray:
        """计算所有原子的parr函数[n,5](N+1,N,N-1,pk-,pk+)"""
        mols=[molN,molP]
        cals=[spin.Calculator(mol) for mol in mols]
        natm=molN.atoms.natm
        vals=np.zeros(shape=(natm,5))
        for c,cal in enumerate(cals):
            vals[:,c]=cal.spin(ctype)
        vals[:,3]=vals[:,0]-vals[:,1]
        vals[:,4]=vals[:,1]-vals[:,2]
        return vals
    
    # 双描述符
    def dual(self,molN:Mole,molP:Mole,ctype:str)->np.ndarray:
        """计算凝聚双描述符[n,4](N+1,N,N-1,f2)
        """
        mols=[molN,self.mol,molP]
        cals=[charge.Calculator(mol) for mol in mols]
        natm=self.mol.atoms.natm
        vals=np.zeros(shape=(natm,4)) # 记录所有原子的电荷
        for c,cal in enumerate(cals):
            cal.form='charge'
            vals[:,c]=cal.charge(ctype)
        vals[:,3]=2*vals[:,1]-vals[:,0]-vals[:,2]
        return vals
    
    def mayerTotalValence(self):
        natm=self.mol.atoms.natm
        vals=np.zeros(natm)
        PS=self.mol.PM@self.mol.SM
        OM=PS*PS.T
        diag=np.diag(PS)
        for i,atom in enumerate(self.mol.atoms):
            u,l=atom.obtBorder
            vals[i]=2*np.sum(diag[u:l])-np.sum(OM[u:l,u:l])
        return vals

    def mayerFreeValence(self):
        pass

    # 电子能差
    def engDiff(self,molN:Mole,molP:Mole)->np.ndarray:
        """计算电子能差，中性到两个状态能量变化，变得越小越好

        Args:
            molN (Mole): N+1 电子
            molP (Mole): N-1 电子

        Returns:
            np.ndarray: 电子能差
        """
        from pywfn.atomprop import energy
        mols=[molN,self.mol,molP]
        cals=[energy.Calculator(mol) for mol in mols]
        natm=self.mol.atoms.natm
        vals=np.zeros(shape=(natm,5))
        for c,cal in enumerate(cals):
            res=cal.atmEngs()
            vals[:,c]=res
        vals[:,3]=vals[:,0]-vals[:,1]
        vals[:,4]=vals[:,1]-vals[:,2]
        return vals

    # 化合价
    def valence(self)->np.ndarray:
        """计算原子化合价，与所有原子相邻的Mayer键级之和"""
        natm=self.mol.atoms.natm
        vals=np.zeros(natm)
        PS=self.mol.PM@self.mol.SM
        OM=PS*PS.T
        for ai,atomi in enumerate(self.mol.atoms):
            ui,li=atomi.obtBorder
            for aj,atomj in enumerate(self.mol.atoms):
                if atomj.idx==atomi.idx:continue
                uj,lj=atomj.obtBorder
                val=np.sum(OM[ui:li,uj:lj]) # 计算键级
                vals[ai]+=val
                if config.SHOW_LEVEL>=1:
                    print(f'{atomi.idx:>3}{atomj.idx:>3}{vals[ai]:>10.4f}')
        return vals

    # 自由价(方向化合价)
    def freeValence(self,atm:int,dirs:np.ndarray)->np.ndarray:
        """计算指定原子的自由价，计算的分子可指定，该原子在不同方向的自由价

        Args:
            atm (int): 原子索引
            mol (Mole | None): 可指定的分子

        Returns:
            np.ndarray: 自由价[d](val)
        """
        from pywfn.bondprop import order
        caler=order.Calculator(self.mol)
        STAND=3.0
        result=[]
        bonds,orders=caler.boundMayer(atm,dirs) # 计算束缚键级
        # print(bmays)
        nebs=self.mol.atom(atm).neighbors # 原子的邻接原子
        for i in range(0,len(orders),len(nebs)):
            valence=STAND-orders[i:i+len(nebs)].sum()
            result.append(valence)
        result=np.array(result)
        return result

    # 亲核亲电自由价 v1
    def neFreeValence_v1(self,atm:int,dirs:np.ndarray,molN:Mole,molP:Mole):
        """计算自由价之差"""
        mols=[molN,self.mol,molP]
        cals=[Calculator(mol) for mol in mols]
        valn,val0,valp=[cal.freeValence(atm,dirs)[1] for cal in cals]
        dirs=val0[:,:-1]
        valsn=(valn[:,-1]-val0[:,-1]).reshape(-1,1)
        valsp=(valp[:,-1]-val0[:,-1]).reshape(-1,1)
        result=np.concatenate([dirs,valsn,valsp],axis=1)
        return result
    
    # 亲核亲电自由价 v2
    def neFreeValence_v2(self,atm:int,dirs:np.ndarray,molN:Mole,molP:Mole):
        """计算自由价之差"""
        mols=[molN,self.mol,molP]
        cals=[Calculator(mol) for mol in mols]
        valn,val0,valp=[cal.freeValence(atm,dirs)[1] for cal in cals]
        val:float=self.valence()[atm-1] # 原子化合价
        dirs=val0[:,:-1]
        valsn=(4-val+(valn[:,-1])-val0[:,-1]).reshape(-1,1)
        valsp=(4-val+(val0[:,-1])-valp[:,-1]).reshape(-1,1)
        result=np.concatenate([dirs,valsn,valsp],axis=1)
        return result

    # 方向福井函数
    def dirFukui(self,atms:list[int],dirs:np.ndarray,molN:Mole,molP:Mole,ctype:str)->np.ndarray:
        """计算方向福井函数

        Args:
            atms (list[int]): 要计算的原子
            dirs (list[np.ndarray]): 原子的方向
            molN (Mole): 多一个电子的分子
            molP (Mole): 少一个电子的分子
            ctype(str): 电荷类型,默认为mulliken

        Returns:
            np.ndarray: 指定方向的福井函数[n,6](N-1,N,N+1,f-,f+,f0)
        """
        
        mols=[molN,self.mol,molP]
        crgs=[mol.multi[0] for mol in mols]
        assert crgs[0]<crgs[1]<crgs[2],"电荷顺序不符"
        calers=[charge.Calculator(mol) for mol in mols]
        
        vals:list[np.ndarray]=[]
        for caler in calers:
            caler.form='number'
            res=caler.dirElectron(atms,dirs,ctype) # 计算方向电子
            vals.append(res)

        natm=self.mol.atoms.num
        result=np.zeros((natm,6))
        result[:,0]=vals[0]
        result[:,1]=vals[1]
        result[:,2]=vals[2]
        result[:,3]=vals[1]-vals[2]
        result[:,4]=vals[0]-vals[1]
        result[:,5]=(vals[0]-vals[2])/2
        return result
    
    # pi双描述符
    def dual_pi(self,molN:Mole,molP:Mole,ctype:str): # 基于π电子的双描述符
        from pywfn.atomprop import charge
        calerN=charge.Calculator(molN)
        calerP=charge.Calculator(molP)
        caler0=charge.Calculator(self.mol)
        calerN.form='charge'
        calerP.form='charge'
        caler0.form='charge'
        piN=calerN.piElectron(ctype)
        piP=calerP.piElectron(ctype)
        pi0=caler0.piElectron(ctype)
        dual=piN[:,-1]+piP[:,-1]-pi0[:,-1]*2
        return dual