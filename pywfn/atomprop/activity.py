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
    def __init__(self,mol:Mol) -> None:
        self.mol=mol

    # 福井函数
    def fukui(self,molN:Mol,molP:Mol,ctype:str)->np.ndarray:
        """计算所有原子的福井函数

        Args:
            molN (Mol): 负电荷分子
            molP (Mol): 正电荷分子
            ctype (Chrgs, optional): 电荷类型. Defaults to 'mulliken'.

        Returns:
            np.ndarray: 福井函数[n,6](N-1,N,N+1,f-,f+,f0)
        """
        mols=[molN,self.mol,molP]
        cals=[charge.Calculator(mol) for mol in mols]
        natm=self.mol.atoms.natm
        vals=np.zeros(shape=(natm,6)) # 记录所有原子的电荷
        for c,cal in enumerate(cals):
            cal.numForm=False
            vals[:,c]=cal.charge(ctype)
        vals[:,3]=vals[:,2]-vals[:,1]
        vals[:,4]=vals[:,1]-vals[:,0]
        vals[:,5]=(vals[:,2]-vals[:,0])/2
        return vals
    
    # parr函数
    def parr(self,molN:Mol,molP:Mol,ctype:str)->np.ndarray:
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
    def dual(self,molN:Mol,molP:Mol,ctype:str)->np.ndarray:
        """计算凝聚双描述符[n,4](N+1,N,N-1,f2)
        """
        mols=[molN,self.mol,molP]
        cals=[charge.Calculator(mol) for mol in mols]
        natm=self.mol.atoms.natm
        vals=np.zeros(shape=(natm,4)) # 记录所有原子的电荷
        for c,cal in enumerate(cals):
            cal.numForm=False
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
    def engDiff(self,molN:Mol,molP:Mol)->np.ndarray:
        """计算电子能差，中性到两个状态能量变化，变得越小越好

        Args:
            molN (Mol): N+1 电子
            molP (Mol): N-1 电子

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
                vals[ai]+=np.sum(OM[ui:li,uj:lj])
        return vals

    # 自由价(方向化合价)
    def freeValence(self,atm:int,dirs:np.ndarray)->np.ndarray:
        """计算指定原子的自由价，计算的分子可指定，该原子在不同方向的自由价

        Args:
            atm (int): 原子索引
            mol (Mol | None): 可指定的分子

        Returns:
            np.ndarray: 自由价[d,5](atm,x,y,z,val)
        """
        from pywfn.bondprop import order
        caler=order.Calculator(self.mol)
        STAND=3.0
        result=[]
        bmays=caler.boundMayer(atm,dirs) # 计算束缚键级
        # print(bmays)
        nebs=self.mol.atom(atm).neighbors # 原子的邻接原子
        for i in range(0,len(bmays),len(nebs)):
            orders=bmays[i:i+len(nebs),-1]
            x,y,z=bmays[i,2:5]
            valence=STAND-sum(orders)
            result.append([atm,x,y,z,valence])
        result=np.array(result)
        return result

    # 亲核亲电自由价 v1
    def neFreeValence_v1(self,atm:int,dirs:np.ndarray,molN:Mol,molP:Mol):
        """计算自由价之差"""
        mols=[molN,self.mol,molP]
        cals=[Calculator(mol) for mol in mols]
        valn,val0,valp=[cal.freeValence(atm,dirs) for cal in cals]
        dirs=val0[:,:-1]
        valsn=(valn[:,-1]-val0[:,-1]).reshape(-1,1)
        valsp=(valp[:,-1]-val0[:,-1]).reshape(-1,1)
        result=np.concatenate([dirs,valsn,valsp],axis=1)
        return result
    
    # 亲核亲电自由价 v2
    def neFreeValence_v2(self,atm:int,dirs:np.ndarray,molN:Mol,molP:Mol):
        """计算自由价之差"""
        mols=[molN,self.mol,molP]
        cals=[Calculator(mol) for mol in mols]
        valn,val0,valp=[cal.freeValence(atm,dirs) for cal in cals]
        val:float=self.valence()[atm-1] # 原子化合价
        dirs=val0[:,:-1]
        valsn=(4-val+(valn[:,-1])-val0[:,-1]).reshape(-1,1)
        valsp=(4-val+(val0[:,-1])-valp[:,-1]).reshape(-1,1)
        result=np.concatenate([dirs,valsn,valsp],axis=1)
        return result

    # 方向福井函数
    def dirFukui(self,atms:list[int],dirs:np.ndarray,molN:Mol,molP:Mol,ctype:str)->np.ndarray:
        """计算方向福井函数

        Args:
            atms (list[int]): 要计算的原子
            dirs (list[np.ndarray]): 原子的方向
            molN (Mol): 多一个电子的分子
            molP (Mol): 少一个电子的分子
            ctype(str): 电荷类型,默认为mulliken

        Returns:
            np.ndarray: 指定方向的福井函数[n,6](N-1,N,N+1,f-,f+,f0)
        """
        mols=[molN,self.mol,molP]
        crgs=[mol.charge for mol in mols]
        # assert crgs[0]<crgs[1]<crgs[2],"电荷顺序不符"
        cals=[charge.Calculator(mol) for mol in mols]
        
        vals:list[np.ndarray]=[]
        for cal in cals:
            cal.numForm=True
            res=cal.dirElectron(atms,dirs,ctype) # 计算方向电子
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
    def dual_pi(self,molN:Mol,molP:Mol,ctype:str): # 基于π电子的双描述符
        from pywfn.atomprop import charge
        calerN=charge.Calculator(molN)
        calerP=charge.Calculator(molP)
        caler0=charge.Calculator(self.mol)
        calerN.numForm=True
        calerP.numForm=True
        caler0.numForm=True
        piN=calerN.piElectron(ctype)
        piP=calerP.piElectron(ctype)
        pi0=caler0.piElectron(ctype)
        dual=piN[:,-1]+piP[:,-1]-pi0[:,-1]*2
        return dual

    
    def onShell(self,shell:Shell):
        from pywfn.utils import printer
        while True:
            printer.options('原子活性',{
                '1':'福井函数',
                '2':'parr函数',
                '3':'双描述符',
                '4':'原子能差',
                '5':'化合价',
                '6':'自由价',
                '7':'方向福井函数',
            })
            opt=input('选择计算活性类型:')
            ctypes={'':'hirshfeld','1':'mulliken','2':'lowdin','3':'hirshfeld','4':'pi_pocv'}
            chrgTip='1.mulliken; 2.lowdin; 3.hirshfeld*; 4.pi_pocv'
            match opt:
                case '1': # 福井函数
                    if len(shell.paths)<3:
                        print('需要至少3个分子')
                        break
                    molN,molP=shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    copt=input('请输入电荷类型: ')
                    if copt not in ctypes.keys():continue
                    result=self.fukui(molN,molP,ctypes[copt])
                    print(f'idx:{"q(N+1)":>10}{"q(N)":>10}{"q(N-1)":>10}{"f-":>10}{"f+":>10}{"f0":>10}')
                    for i,(en,e0,ep,fn,fp,f0) in enumerate(result):
                        print(f'{i+1:>3d}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{fn:>10.4f}{fp:>10.4f}{f0:>10.4f}')
                case '2': # parr函数
                    molN,molP=shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    opt=input('请输入电荷类型:')
                    if opt not in ctypes.keys():continue
                    result=self.parr(molN,molP,ctypes[opt])
                    print(f'idx:{"s(N+1)":>10}{"s(N-1)":>10}{"f-":>10}{"f+":>10}{"f0":>10}')
                    for i,(sn,s0,sp,ve,vn) in enumerate(result):
                        print(f'{i+1:>3d}:{sn:>10.4f}{s0:>10.4f}{sp:>10.4f}{ve:>10.4f}{vn:>10.4f}')
                case '3': # 双描述符
                    if len(shell.paths)<3:
                        print('需要至少3个分子')
                        break
                    molN,molP=shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    opt=input('请输入电荷类型: ')
                    if opt not in ctypes.keys():continue
                    result=self.dual(molN,molP,ctypes[opt])
                    print(f'idx:{"q(N+1)":>10}{"q(N)":>10}{"q(N-1)":>10}{"CDD":>10}')
                    for i,(en,e0,ep,val) in enumerate(result):
                        print(f'{i+1:>3}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{val:>10.4f}')
                case '4': # 原子能差
                    molN,molP=shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    result=self.engDiff(molN,molP)
                    print(f'idx:{"E(N+1)":>10}{"E(N)":>10}{"E(N-1)":>10}{"d-":>10}{"d+":>10}')
                    for i,(en,e0,ep,ve,vn) in enumerate(result):
                        print(f'{i+1:>3d}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{ve:>10.4f}{vn:>10.4f}')
                case '5': # 化合价
                    result=self.valence()
                    for i,val in enumerate(result):
                        print(f'{i+1:>3d}:{val:>10.4f}')
                case '6': # 自由价，可以输入方向或内置方向
                    from pywfn.atomprop import direction
                    dirCaler=direction.Calculator(self.mol)
                    atm=shell.input.Integ(tip='输入原子编号: ',count=1)[0]
                    dirs=shell.input.Float(tip='?指定投影方向: ',count=3)
                    if dirs is None:
                        dirs=dirCaler.reactions(atm)
                    else:
                        dirs=np.array(dirs).reshape(1,3)
                    result=self.freeValence(atm,dirs)
                    for a,x,y,z,v in result:
                        print(f'{atm:>3d}:{x:>10.4f}{y:>10.4f}{z:>10.4f}{v:>10.4f}')
                case '7': # 方向fukui函数
                    # self.mols=shell.input.Moles(num=3)
                    copt=input('请输入电荷类型: ')
                    if copt not in ctypes.keys():continue
                    atms=shell.input.Integ(tip='输入原子编号: ')
                    assert atms is not None,"输入错误"
                    dirs=shell.input.Float(tip='?指定投影方向: ',count=3)
                    dirs=np.array(dirs).reshape(1,3)
                    molN,molP=shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    result=self.dirFukui(atms,dirs,molN,molP,ctypes[copt])
                    for i,val in enumerate(result):
                        print(f'{i+1:>3d} {val:>8.4f}')
                case _:
                    break