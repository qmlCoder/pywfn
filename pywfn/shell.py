"""
模仿multiwfn以命令窗口的形式操作程序计算各种想要的东西
可以彩色输出
绿色-结果
黄色-警告
红色-错误
在没有特别提示的地方,空输入代表返回上级
功能模块划分与pywfn设计一致
"""

from pathlib import Path
import re, sys

from pywfn.reader import get_reader
from pywfn.utils import printer
from pywfn import data
from pywfn.base import Mole
from pywfn import config
import numpy as np


class Shell:
    def __init__(self):
        self.paths: list[str] = []
        self.mols: dict[str,Mole] = {} # 路径:分子的对应
        config.IF_SHELL=True
        printer.info(data.start)
        self.input = Inputer(self)
        self.atomPage=AtomPage(self)
        self.bondPage=BondPage(self)
        self.molePage=MolePage(self)
        self.gridPage=GridPage(self)
        self.toolPage=ToolPage(self)
        self.savePage=SavePage(self)
        self.chemPage=ChemPage(self)

    def homePage(self):
        opts = {
            "0": "导入文件",
            "1": "原子性质",
            "2": "键の性质",
            "3": "分子性质",
            "4": "空间性质",
            "5": "实用工具",
            "6": "文件转换",
            "7": "化学函数"
        }
        while True:
            opts["0"]=f"导入文件({len(self.paths)})"
            printer.options('首页', opts)
            opt = input("输入功能选项: ")
            match opt:
                case "0":
                    self.input.Files()
                case "1":
                    self.atomPage.home()
                case "2":
                    self.bondPage.home()
                case "3":
                    self.molePage.home()
                case "4":
                    self.gridPage.home()
                case "5":
                    self.toolPage.home()
                case "6":
                    self.savePage.home()
                case "7":
                    self.chemPage.home()
                case 'q':
                    break
                case _:
                    printer.warn("命令不存在")



class AtomPage:
    def __init__(self,shell:Shell) -> None:
        self.name='原子属性'
        self.opts={
            '1':'电子布局',
            '2':'活性指标',
            '3':'原子能量',
            
        }
        self.shell=shell
    
    def home(self):
        while True:
            printer.options(self.name,self.opts)
            opt=input('请输入对应选项：')
            match opt:
                case '1':
                    self.charge()
                case '2':
                    self.activity()
                case '3':
                    self.energy()
                case '4':
                    from pywfn.atomprop import activity
                    mol=self.shell.input.Moles(tip='输入要计算活性的分子',num=1)[0]
                    caler=activity.Calculator(mol)
                case _:
                    break
    
    def charge(self):
        from pywfn.atomprop import charge
        caler=charge.Calculator(self.shell.input.Moles()[0])
        chrgMap={'':'mulliken','1':'mulliken','2':'lowdin','3':'space','4':'hirshfeld'}
        chrgStr='1. Mulliken[*]; 2. lowdin; 3. sapce; 4. hirshfeld'
        while True:
            printer.options('原子电荷',{
                '1':'mulliken电荷',
                '2':'lowdin电荷',
                '3':'空间积分电荷',
                '4':'hirshfeld电荷',
                '5':'方向电子分布',
                '6':'π 电子分布',
                '7':'电子自旋分布'
            })
            opt=input('请选择计算类型: ')
            match opt:
                case '1': # Mulliken
                    charges=caler.mulliken()
                    for i,val in enumerate(charges):
                        print(f'{i+1:>3d}: {val:>8.4f}')
                    print(f'sum:{np.sum(charges)}')

                case '2': # Lowdin
                    charges=caler.lowdin()
                    for i,val in enumerate(charges):
                        print(f'{i+1:>3d}: {val:>8.4f}')
                    print(f'sum:{np.sum(charges)}')
                
                case '3': # 空间积分电荷
                    charges=caler.sapce()
                    for i,val in enumerate(charges):
                        print(f'{i+1:>3d}: {val:>8.4f}')
                    print(f'sum:{np.sum(charges)}')

                case '4': # hirshfeld电荷
                    charges=caler.hirshfeld()
                    for i,val in enumerate(charges):
                        print(f'{i+1:>3d}: {val:>8.4f}')
                    print(f'sum:{np.sum(charges)}')

                case '5': # 方向电子
                    print(chrgStr)
                    opt=input('选择电荷类型: ')
                    if opt not in chrgMap.keys():return
                    chrg=chrgMap[opt]
                    # numStr=input('输入原子编号: ')
                    atm=self.shell.input.Integ(tip='输入原子编号: ',count=1)[0]
                    dir=self.shell.input.Float(tip='输入原子向量: ',count=3)
                    assert dir is not None,'原子向量输入不正确'
                    dirs=np.array(dir).reshape(1,3)
                    elects=caler.dirElectron(atms=[atm],dirs=dirs,ctype=chrg)
                    ele=elects[atm-1]
                    x,y,z=dir
                    print(f'{atm:>3d} ({x:>6.2f},{y:>6.2f},{z:>6.2f}):{ele:>8.4f}')
                
                case '6': # pi电子
                    print(chrgStr)
                    opt=input('选择电荷类型: ')
                    if opt not in chrgMap.keys():return
                    chrg=chrgMap[opt]
                    elects=caler.piElectron(chrg)
                    for i,(idx,x,y,z,val) in enumerate(elects):
                        print(f'{i+1:>3d}:{val:>8.4f}')
                    print(f'sum:{elects[:,-1].sum()}')
                
                case '7':
                    print(chrgStr)
                    opt=input('选择电荷类型: ')
                    spins=caler.spin(chrg='mulliken')
                    for i,spin in enumerate(spins):
                        print(f'{i+1:>3d}:{spin:>10.4f}')
                    print(f'自旋之和: {sum(spins):>10.4f}')

                case _:
                    break
    
    def activity(self):
        from pywfn.atomprop import activity
        mol=self.shell.input.Moles(tip='输入要计算活性的分子',num=1)[0]
        caler=activity.Calculator(mol)
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
                    if len(self.shell.paths)<3:
                        print('需要至少3个分子')
                        break
                    molN,molP=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    copt=input('请输入电荷类型: ')
                    if copt not in ctypes.keys():continue
                    result=caler.fukui(molN,molP,ctypes[copt])
                    print(f'idx:{"q(N+1)":>10}{"q(N)":>10}{"q(N-1)":>10}{"f-":>10}{"f+":>10}{"f0":>10}')
                    for i,(en,e0,ep,fn,fp,f0) in enumerate(result):
                        print(f'{i+1:>3d}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{fn:>10.4f}{fp:>10.4f}{f0:>10.4f}')
                case '2': # parr函数
                    molN,molP=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    opt=input('请输入电荷类型:')
                    if opt not in ctypes.keys():continue
                    result=caler.parr(molN,molP,ctypes[opt])
                    print(f'idx:{"s(N+1)":>10}{"s(N-1)":>10}{"f-":>10}{"f+":>10}{"f0":>10}')
                    for i,(sn,s0,sp,ve,vn) in enumerate(result):
                        print(f'{i+1:>3d}:{sn:>10.4f}{s0:>10.4f}{sp:>10.4f}{ve:>10.4f}{vn:>10.4f}')
                case '3': # 双描述符
                    if len(self.shell.paths)<3:
                        print('需要至少3个分子')
                        break
                    molN,molP=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    print(chrgTip)
                    opt=input('请输入电荷类型: ')
                    if opt not in ctypes.keys():continue
                    result=caler.dual(molN,molP,ctypes[opt])
                    print(f'idx:{"q(N+1)":>10}{"q(N)":>10}{"q(N-1)":>10}{"CDD":>10}')
                    for i,(en,e0,ep,val) in enumerate(result):
                        print(f'{i+1:>3}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{val:>10.4f}')
                case '4': # 原子能差
                    molN,molP=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    result=caler.engDiff(molN,molP)
                    print(f'idx:{"E(N+1)":>10}{"E(N)":>10}{"E(N-1)":>10}{"d-":>10}{"d+":>10}')
                    for i,(en,e0,ep,ve,vn) in enumerate(result):
                        print(f'{i+1:>3d}:{en:>10.4f}{e0:>10.4f}{ep:>10.4f}{ve:>10.4f}{vn:>10.4f}')
                case '5': # 化合价
                    result=caler.valence()
                    for i,val in enumerate(result):
                        print(f'{i+1:>3d}:{val:>10.4f}')
                case '6': # 自由价，可以输入方向或内置方向
                    from pywfn.atomprop import direction
                    mol=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=1)[0]
                    dirCaler=direction.Calculator(mol)
                    atm=self.shell.input.Integ(tip='输入原子编号: ',count=1)[0]
                    dirs=self.shell.input.Float(tip='?指定投影方向: ',count=3)
                    if dirs is None:
                        dirs=dirCaler.reactions(atm)
                    else:
                        dirs=np.array(dirs).reshape(1,3)
                    result=caler.freeValence(atm,dirs)
                    for a,x,y,z,v in result:
                        print(f'{atm:>3d}:{x:>10.4f}{y:>10.4f}{z:>10.4f}{v:>10.4f}')
                case '7': # 方向fukui函数
                    # self.mols=shell.input.Moles(num=3)
                    copt=input('请输入电荷类型: ')
                    if copt not in ctypes.keys():continue
                    atms=self.shell.input.Integ(tip='输入原子编号: ')
                    assert atms is not None,"输入错误"
                    dirs=self.shell.input.Float(tip='?指定投影方向: ',count=3)
                    dirs=np.array(dirs).reshape(1,3)
                    molN,molP=self.shell.input.Moles(tip='分别输入N+1和N-1个电子的分子',num=2)
                    result=caler.dirFukui(atms,dirs,molN,molP,ctypes[copt])
                    for i,val in enumerate(result):
                        print(f'{i+1:>3d} {val:>8.4f}')
                case _:
                    break

    def energy(self):
        from pywfn.atomprop import energy
        caler=energy.Calculator(self.shell.input.Moles()[0])
        while True:
            printer.options('原子能量',{
                '1':'原子  电子能 (将分子轨道能量分配到原子上)',
                '2':'原子pi电子能 (pi分子轨道能量分配到原子上)'
            })
            opt=input('输入原子能类型:')
            if opt=='1':
                energy=caler.atmEngs()
                for i,eng in enumerate(energy):
                    print(f'{i+1:>3}:{eng:>10.4f}')
            elif opt=='2':
                engs=caler.atmPiEngs()
                for i,eng in enumerate(engs):
                    print(f'{i+1:>3}:{eng:>10.4f}')
            else:
                break


class BondPage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        from pywfn.bondprop import order
        printer.options('键の属性',{
        '1':'各种键级'
        })
        opt=input('请输入要计算的键属性:')
        if opt=='1':
            self.order()
            
    
    def order(self):
        from pywfn.bondprop import order
        mols=self.shell.input.Moles()
        caler=order.Calculator(mols[0])
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
                    orders=caler.mayer()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '2':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        
                        a1,a2=opt.split('-')
                        dirs=self.shell.input.Float(tip='输入方向: ',count=3)
                        dirs=np.array(dirs).reshape(1,3)
                        result=caler.dirMayer(bond=[int(a1),int(a2)],dirs=dirs)
                        for a1,a2,x,y,z,val in result:
                            print(f'{int(a1):>2d}-{int(a2):>2d}({x:>8.4f} {y:>8.4f} {z:>8.4f}):{val:>8.4f}')
                case '3':
                    orders=caler.pi_pocv()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '4':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        a1,a2=opt.split('-')
                        order=caler.pi_smo(bond=[int(a1),int(a2)])
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{order:>8.4f}')
                case '5':
                    orders=caler.hmo()
                    for a1,a2,val in orders:
                        print(f'{int(a1):>2d}-{int(a2):>2d}:{val:>8.4f}')
                case '6':
                    while True:
                        opt=input('请输入需要计算的键，例如(1-2): ')
                        if not opt:break
                        a1,a2=opt.split('-')
                        orders=caler.decompose([int(a1),int(a2)])
                        sig,piz,pix,det=orders
                        print(f'σ : {sig:>10.4f}')
                        print(f'πz: {piz:>10.4f}')
                        print(f'πx: {pix:>10.4f}')
                        print(f'δ : {det:>10.4f}')
                case _:
                    break
        

class MolePage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        printer.options('分子属性',{
            '1':'芳香性',
        })
        opt=input('请输入对应选项：')
        match opt:
            case '1':
                self.aromacity()

    def aromacity(self):
        from pywfn.utils import printer
        from pywfn.moleprop import aromatic
        mol=self.shell.input.Moles()[0]
        caler=aromatic.Calculator(mol)
        while True:
            printer.options('芳香性',{
                '1':'piSD  (根据pi键级的标准差)',
                '2':'piMSD (根据pi键级的均值和标准差)',
                '3':'piMED',
                '4':'HOMED',
            })
            opt=input('输入芳香性类型:')
            
            match opt:
                case '1':
                    ring=self.shell.input.Integ('?输入环编号: ')
                    if len(ring)==0:ring=None
                    result=caler.pisd(ring)
                    print(f'{result}')
                case '2':
                    result=caler.pimsd(ring=ring,ratio=0.5)
                    print(f'{result}')
                case '3':
                    result=caler.pimed()
                    print(f'{result}')
                case '4':
                    result=caler.homed()
                    print(f'{result}')
                case _:
                    break


class GridPage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        from pywfn.gridprop import read_grid
        while True:
            printer.options('空间性质',{
                '1':'  波函数',
                '2':'电子密度',
                '3':'  静电势',
            })
            ctype=input('请选择计算的性质: ')
            mol=self.shell.input.Moles()[0]
            show_dict({'1':'直线采点','2':'平面采点','3':'空间采点','4':'分子空间','5':'地图映射','6':'范德华表面'})
            gidx=input('请选择采点方式: ')
            obje=read_grid(self.shell,gidx,mol)
            match ctype:
                case '1':
                    self.wfnfunc(mol,obje)
                case '2':
                    self.density(mol,obje)
                case '3':
                    self.potential(mol,obje)
    
    def wfnfunc(self,mol,obje):
        from pywfn.gridprop import wfnfunc,save_file
        caler=wfnfunc.Calculator(mol)
        atms=self.shell.input.Integ(tip='?输入要计算的原子: ')
        if not atms:atms=mol.atoms.atms
        obts=self.shell.input.Integ(tip='*输入要计算的轨道: ')
        wfns=caler.atmWfns(obje.grid,atms,obts)
        save_file(caler,obje,wfns,obts=obts)
    
    def density(self,mol,obje):
        from pywfn.gridprop import density,save_file
        caler=density.Calculator(mol)
        atms=self.shell.input.Integ(tip='?输入要计算的原子: ')
        if not atms:atms=mol.atoms.atms
        dens=caler.atmDens(obje.grid,atms)
        dens=dens.reshape(1,-1)
        save_file(caler,obje,dens)
    
    def potential(self,mol,obje):
        from pywfn.gridprop import potential,save_file
        caler=potential.Calculator(mol)
        pots=caler.molPotential(obje.grid)
        save_file(caler,obje,pots)


class ToolPage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        from pathlib import Path
        cwd=Path.cwd()
        printer.options('实用工具',{
            '1':'提取  SI 信息',
            '2':'分割SCAN 文件',
            '3':'分割 IRC 文件',
            '4':'分割link 任务',
            '5':'拼接 gjf 文件',
            '6':'gjf 环心添加Bq原子',
            '7':'gjf 分子增删电子数'
        })
        opt=input('请输入选项：')
        match opt:
            case '1': # 提取SI信息
                from pywfn.tools import extractSI
                paths=self.shell.input.Paths()
                tool=extractSI.Tool(paths)
                path0=Path(paths[0]) # 第一个分子
                spath=path0.parent/'SI.txt' # 第一个文件所在的文件夹
                tool.save(f'{spath}')
            case '2': # 分割SCAN文件
                from pywfn.tools import spiltScan
                paths=self.shell.input.Paths()
                for path in paths:
                    spiltScan.Tool(path).save()
            case '3': # 分割IRC文件
                from pywfn.tools import splitIrc
                paths=self.shell.input.Paths()
                for path in paths:
                    splitIrc.Tool(path).split()
            case '4': # 分割link任务
                from pywfn.tools import splitLink
                paths=self.shell.input.Paths()
                for path in paths:
                    splitLink.Tool(path).split()
            case '5': # 拼接gjf文件
                from pywfn.tools import joinGjf
                paths=self.shell.input.Paths()
                joinGjf.Tool(paths).save(f'{cwd}/join.gjf')
            case '6': # 环心添加Bq原子
                from pywfn.tools import editGjf
                printer.info('在gjf文件指定环的中心添加Bq原子，方便NICS计算')
                mol=self.shell.input.Moles(num=1)[0]
                tool=editGjf.Tool(mol)
                rings=input('输入环编号: ')
                rings=[[int(atm) for atm in ring] for ring in rings.split(';')]
                tool.addRingBq(rings)
                path=Path(mol.reader.path)
                path=(path.parent/f'{path.stem}_ringBq.gjf')
                tool.save(f'{path}')
            case '7': # 加/减gjf电子数
                from pywfn.tools import editGjf
                mols=self.shell.input.Moles()
                nele=input('输入加减电子数: ')
                nele=int(nele)
                if nele>0:
                    sufx='n'
                elif nele<0:
                    sufx='p'
                else:
                    raise ValueError('加减电子数不能为0')
                for mol in mols:
                    tool=editGjf.Tool(mol)
                    tool.addElectron(nele)
                    path=Path(mol.reader.path)
                    tool.save(f'{path.parent}/{path.stem}_{sufx}{abs(nele)}.gjf')


class SavePage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        from pywfn.writer import GjfWriter,MolWriter,XyzWriter
        from pywfn.data import temps
        printer.options('文件转换',{
            '0':'导出模板',
            '1':'gjf文件',
            '2':'mol文件',
            '3':'xyz文件',
        })
        opt=input('请选择文件类型:')
        
        match opt:
            case '0':
                (Path().cwd()/'temp.gjf.txt').write_text(temps.gjf)
                (Path().cwd()/'temp.mol.txt').write_text(temps.mol)
                (Path().cwd()/'temp.mol.xyz').write_text(temps.xyz)
            case '1':
                mols=self.shell.input.Moles()
                for mol in mols:
                    writer=GjfWriter().fromMol(mol)
                    path=Path(mol.reader.path).with_suffix('.gjf')
                    writer.save(f'{path}')
            case '2':
                mols=self.shell.input.Moles()
                for mol in mols:
                    writer=MolWriter().fromMol(mol)
                    path=Path(mol.reader.path).with_suffix('.mol')
                    writer.save(f'{path}')
            case '3':
                mols=self.shell.input.Moles()
                for mol in mols:
                    writer=XyzWriter().fromMol(mol)
                    path=Path(mol.reader.path).with_suffix('.xyz')
                    writer.save(f'{path}')
            case _:
                return


class ChemPage:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell
    
    def home(self):
        from pywfn import chems
        printer.options('化学计算器',{
            '1':'计算反应活化能',
            '2':'计算立体选择ee值',
            '3':'计算玻尔兹曼分布',
        })
        opt=input('请输入对应选项：')
        match opt:
            case '1':
                K=float(input('催化常数K:'))
                T=float(input('反应温度T:'))
                print('反应活化能:',chems.reaActEne(K,T),'kJ/mol')
            case '2':
                deR=float(input('反应热力学参数:'))
                deS=float(input('立体选择热力学参数:'))
                T=float(input('反应温度T:'))
                print('立体选择ee值:',chems.steSelEE(deR,deS,T))
            case '3':
                engs=input('能量列表(空格分隔):').split()
                engs=[float(i) for i in engs]
                engs=np.array(engs)
                print('玻尔兹曼分布:',chems.bezm(engs))
            case _:
                return


class Inputer:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell

    def input(self, tip: str):
        res = input(tip)
        if res == "q":
            sys.exit()
        return res
    
    def Float(self,tip:str,count:int=0)->list[float]|None:
        """输入浮点数"""
        while True:
            numStr=input(tip)
            if numStr=="":return None
            nums=[float(e) for e in numStr.split(',')]
            if len(nums)!=count:continue
            return nums

    def Integ(self,tip:str,count:int=0)->list[int]:
        """输入整数"""
        def intSplit(frg: str):
            if "-" in frg:
                num1, num2 = [int(e) for e in frg.split("-")]
                assert num1 < num2, "数字1要大于数字2"
                return list(f"{e}" for e in range(num1, num2 + 1))
            else:
                return [frg]

        while True:
            numStr=input(tip)
            if numStr=="":
                return []
            frgs = numStr.split(',')
            nums = []
            for frg in frgs:
                nums += intSplit(frg)
            nums = [int(e) for e in nums]
            if count!=0 and len(nums)!=count:
                printer.warn(f"长度不符合要求,应该为{count}")
                continue
            return nums


    def Vector(self, tip: str, norm: bool = False, must=False) -> np.ndarray|None:
        """
        输入向量
        tip:提示
        norm:是否归一化
        must:是否必须
        """
        while True:
            opt = input(tip)
            if opt == "":
                if must:
                    continue
                else:
                    return None
            nums = re.split(r"[,，]", opt)
            if len(nums) != 3:
                continue
            nums = np.array(nums, dtype=np.float32)
            if norm:
                nums /= np.linalg.norm(nums)
            return nums

    @staticmethod
    def Bool(tip: str, default: bool = True):
        opt = input(tip)
        match opt:
            case "":
                return default
            case "y":
                return True
            case "n":
                return False
            case _:
                printer.warn("请输入正确选项!")
                return Inputer.Bool(tip, default)

    def Files(self):
        """
        用户输入文件
        """
        path=input('请输入文件路径: ')
        types = [".log", ".out", ".fch", ".gjf"]  # 支持的文件类型
        pathObj = Path(path)
        if pathObj.is_file():  # 如果是文件
            if pathObj.suffix in types:
                printer.info(f"共1个文件")
                self.shell.paths.append(path)
            else:
                printer.warn("不支持的文件类型")
        elif pathObj.is_dir():  # 如果是文件夹
            count=0
            for each in pathObj.iterdir():  # 对文件夹中的每个文件进行循环
                if each.suffix in types:  # 如果文件类型是支持的文件类型
                    self.shell.paths.append(f'{each}')
                    count+=1
            printer.info(f"输入{count}个文件")

    def Paths(self,count:int=0)->list[str]:
        """
        用户选择文件
        count控制数量,如果count为0则没有限制
        """
        paths = self.shell.paths
        printer.info(f"共{len(paths)}个文件:")
        printer.bar()
        for p,path in enumerate(paths):
            print(f'{p:>2} {path}')

        while True:
            if len(paths)==0:
                printer.warn("你还没有输入任何文件!!!")
                self.Files()
                continue
            if len(paths)<count:
                printer.warn(f"当前文件数量不够,至少需要{count}个")
                self.Files()
                continue
            idxs = self.Integ(f"输入分子编号: ")
            if len(idxs) != count and count!=0:
                printer.warn(f"需要输入{count}个数字,当前{len(idxs)}个")
            else:
                break
            if count==0:break
        return [paths[int(idx)] for idx in idxs]

    def Moles(self,tip:str='',num:int=0) -> list[Mole]:
        """
        获取当前文件的分子，列举出当前读取的文件让用户选择，将用户选择的文件对应为分子
        num:分子的数量
        mtype:返回路径还是分子类
        """
        if tip!='':print(f'{tip}')
        if num!=0:print(f'请输入{num}个分子')
        paths=self.Paths(num)
        
        mols = []
        for path in paths:
            if path not in self.shell.mols.keys():
                self.shell.mols[path] = Mole(get_reader(path))
            mols.append(self.shell.mols[path])
        
        return mols
    

def show_dict(dict:dict):
    for k,v in dict.items():
        print(f'{k}: {v}; ',end='')
    print('\n')