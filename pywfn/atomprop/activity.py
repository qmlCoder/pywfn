"""
计算原子的反应活性，福井函数、parr函数等，一般需要多个分子才行
"""
from pywfn.base import Mole
from pywfn.atomprop import charge
import numpy as np
from pywfn import config
from pywfn import core

class Calculator:
    def __init__(self,mole:Mole) -> None:
        self.mole=mole
        self.caler=core.atomprop.activity.Calculator(mole.mole) # type: ignore # 核心计算器

    # 福井函数
    def fukui(self,mole_n:Mole,mole_p:Mole,ctype:str="mulliken"):
        """计算所有原子的福井函数

        Args:
            molN (Mole): 负电荷分子
            molP (Mole): 正电荷分子
            ctype (Chrgs, optional): 电荷类型. Defaults to 'mulliken'.

        Returns:
            np.ndarray: 福井函数[n,7](N,N+1,N-1,f-,f+,f0,df)
        """
        vals=self.caler.fukui(mole_n.mole,mole_p.mole,ctype)
        return np.array(vals)

    # 活性矢量
    def vector(self,atm:int,dir:list[float]):
        """计算指定原子的自由价，计算的分子可指定，该原子在不同方向的自由价

        Args:
            atm (int): 原子索引

        Returns:
            np.ndarray: 自由价[d](val)
        """
        return self.caler.vector(atm,dir)

    # 某个原子的方向福井函数
    def fukui_dir(self,atm:int,dir:list[float],mole_n:Mole,mole_p:Mole,ctype:str="mulliken"):
        """计算方向福井函数

        Args:
            atms (int): 要计算的原子
            dirs (np.ndarray): 原子的方向，可以为多个[n,3]
            molN (Mole): 多一个电子的分子
            molP (Mole): 少一个电子的分子
            ctype(str): 电荷类型,默认为mulliken

        Returns:
            np.ndarray: 指定方向的福井函数[n,6](N,N+1,N-1,f-,f+,f0,df)
        """
        vals=self.caler.fukui_dir(atm,dir,mole_n.mole,mole_p.mole,ctype)
        return np.array(vals)
        
    
    def fukui_pi(self,mole_n:Mole,mole_p:Mole,ctype:str='mulliken'):
        """基于pi电子数的福井函数"""
        dirs,vals=self.caler.fukui_pi(mole_n.mole,mole_p.mole,ctype)
        return dirs, np.array(vals)

    # parr函数
    def parr(self,molN:Mole,molP:Mole,ctype:str)->np.ndarray:
        """计算所有原子的parr函数[n,5](N+1,N,N-1,pk-,pk+)"""
        mols=[molN,molP]
        cals=[charge.Calculator(mol) for mol in mols]
        natm=molN.atoms.natm
        vals=np.zeros(shape=(natm,2))
        for c,cal in enumerate(cals):
            vals[:,c]=cal.spin(ctype)
        return vals
    
    def mayerTotalValence(self):
        natm=self.mole.atoms.natm
        vals=np.zeros(natm)
        PS=self.mole.PM@self.mole.SM
        OM=PS*PS.T
        diag=np.diag(PS)
        for i,atom in enumerate(self.mole.atoms):
            u,l=atom.obtBorder
            vals[i]=2*np.sum(diag[u:l])-np.sum(OM[u:l,u:l])
        return vals

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
        mols=[molN,self.mole,molP]
        cals=[energy.Calculator(mol) for mol in mols]
        natm=self.mole.atoms.natm
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
        natm=self.mole.atoms.natm
        vals=np.zeros(natm)
        PS=self.mole.PM@self.mole.SM
        OM=PS*PS.T
        for ai,atomi in enumerate(self.mole.atoms):
            ui,li=atomi.obtBorder
            for aj,atomj in enumerate(self.mole.atoms):
                if atomj.idx==atomi.idx:continue
                uj,lj=atomj.obtBorder
                val=np.sum(OM[ui:li,uj:lj]) # 计算键级
                vals[ai]+=val
                if config.SHOW_LEVEL>=1:
                    print(f'{atomi.idx:>3}{atomj.idx:>3}{vals[ai]:>10.4f}')
        return vals