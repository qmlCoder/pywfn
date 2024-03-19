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
import re,sys

from pywfn.atomprop import mullikenCharge
from pywfn.reader import get_reader
from pywfn import tools
from pywfn import config
from pywfn.utils import printer
from pywfn import data
from pywfn.base import Mol
import numpy as np

class Shell:
    def __init__(self):
        self.paths:list[Path]=None
        printer.ifShell=True
        printer.info(data.start)
        self.input=Input()
        

    def homePage(self):
        opts=[
            ['0','导入文件'],
            ['1','键级计算'],
            ['2','原子属性'],
            ['3','实用工具'],
            ['4','导出文件'],
            ['5','程序设置']
        ]
        while True:
            opt=self.input.Option('首页',opts,canNone=False)
            match opt:
                case '0':self.paths=self.input.Files()
                case '1':self.bondProp()
                case '2':self.atomProp()
                case '3':self.toolsPage()
                case '4':self.writePage()
                case '5':self.setingPage()
                case _:printer.warn('命令不存在')
        
    def atomProp(self):
        """计算原子属性"""
        opts=[
            ['1','Mulliken 电荷分布'],
            ['2','Mulliken 电子自旋'],
            ['3','pi 电子分布'],
            ['4','方向 电子分布'],
            ['5','方向 电子自旋'],
            ['6','原子 自由价'],
            ['7','原子 相关方向'],
        ]
        path=self.paths[0]
        mol=Mol(reader=get_reader(path))
        while True:
            opt=self.input.Option('计算原子属性',opts)
            match opt:
                case None:break
                
                case '1': # mulliken 电荷分布
                    caler=mullikenCharge.Calculator(mol)
                    caler.print(caler.resStr())
                case '2': # mulliken 电子自旋
                    from pywfn.atomprop import mullikenSpin
                    caler=mullikenSpin.Calculator(mol)
                    caler.print(caler.resStr())
                case '3': # π 电子分布
                    from pywfn.atomprop import piElectron
                    caler=piElectron.Calculator(mol)
                    caler.print(caler.resStr())
                case '4': # 方向 电子分布
                    from pywfn.atomprop import directElectron
                    caler=directElectron.Calculator(mol)
                    while True:
                        atom=self.input.Number(tip='请输入要计算的原子: ',type_='int',length=1)
                        if atom is None:break
                        vect=self.input.Number(tip='请输入要计算的方向: ',length=3)
                        if vect is None:break
                        caler.atoms=[atom]
                        caler.vects=[np.array(vect)]
                        caler.printRes()
                case '5': # 方向 电子自旋
                    from pywfn.atomprop import directSpin
                    caler=directSpin.Calculator(mol)
                    while True:
                        atom=self.input.Number(tip='请输入要计算的原子: ',type_='int',length=1)
                        if atom is None:break
                        vect=self.input.Number(tip='请输入要计算的方向: ',length=3)
                        if vect is None:break
                        caler.atoms=[atom]
                        caler.vects=[np.array(vect)]
                        caler.printRes()
                case '6': # 原子自由价
                    from pywfn.atomprop import freeValence
                    caler=freeValence.Calculator(mol)
                    while True:
                        idx=self.input.Number(tip='请输入原子编号: ',type_='int',length=1)
                        if idx is None:break
                        caler.direct=self.input.Number('输入方向向量(,)[法向量]: ',length=3)
                        caler.print(caler.resStr([idx]))
                case '7':
                    while True:
                        idx=self.input.Number(tip='请输入原子编号: ',type_='int',length=1)
                        if idx is None:break
                        atom=mol.atom(idx)

                        normal=atom.get_Normal()
                        obtWay=atom.get_obtWay(len(mol.O_obts)-1)

                        if normal is not None:
                            printer.vector('原子法向量',normal)
                        else:
                            printer.warn('原子没有法向量')
                        printer.vector('HOMO轨道方向',obtWay)
                        
                        
                case _:
                    return
                         
    def bondProp(self):
        """计算各种键级"""
        opts=[
            ['1','π 键级(OP)'],
            ['2','π 键级(SMO)'],
            ['3','Mayer 键级'],
            ['4','HMO 键级'],
        ]
        if len(self.paths)>1:printer.info('计算第一个分子')
        file=self.paths[0]
        reader=get_reader(file)
        mol=Mol(reader)
        while True:
            match self.input.Option('计算各种键级',opts):
                case None:break
                case '1':
                    from pywfn.bondorder import piDM
                    caler=piDM.Calculator(mol)
                case '2':
                    from pywfn.bondorder import piSM
                    caler=piSM.Calculator(mol)
                case '3':
                    from pywfn.bondorder import mayer
                    caler=mayer.Calculator(mol)
                case '4':
                    from pywfn.bondorder import hmo
                    caler=hmo.Calculator(mol)
                case _:
                    printer.warn('选项不存在')
                    continue
            while True:
                bond=self.input.Number('输入两原子的编号(,)[]: ',type_='int',length=2)
                if bond is None:break
                caler.bond=bond
                caler.print()

    def toolsPage(self):
        """进入实用工具页面"""
        opts=[
            ['1','分割 扫描 文件'],
            ['2','分割 IRC 文件'],
            ['3','分割 link 任务'],
            ['4','生成 PES 文件'],
            ['5','计算 反应活化能'],
            ['6','计算 立体选择性ee值'],
        ]
        while True:
            opt=self.input.Option('使用工具页面',opts)
            match opt:
                case None:break
                case '1':
                    for i,path in enumerate(self.paths):
                        tools.SplitScan(path).split()
                    printer.res('文件分割完成 >_<')
                case '2':
                    for i,path in enumerate(self.paths):
                        tools.SplitIrc(path).split()
                    printer.res('文件分割完成 >_<')
                case '3':
                    for i,path in enumerate(self.paths):
                        tools.SplitLink(path).split()
                case '4':
                    from pywfn.tools.getPES import Tool
                    self.showFiles()
                    tool=Tool(self.paths)
                    tool.route=self.input.Number(tip='路径索引：',type_='int')
                    tool.create()
                case '5':
                    K=self.input.Number(tip='输入K(催化常数): ',length=1)
                    T=self.input.Number(tip='输入T(反应温度,K): ',length=1)
                    R=8.314
                    c=2*1e10*T
                    x=np.log(K/c)
                    g0=-x*R*T*1e-3 # kJ/mol
                    g1=g0/4.184
                    printer.res(f'G(吉布斯自由能):\n{g0:>10.4f} kJ/mol\n{g1:>10.4f} Kcal/mol')
                case '6':
                    deR,deS=self.input.Number(tip='分别输入ΔG(R),ΔG(S): ',length=2)
                    T=self.input.Number(tip='输入T(反应温度,K): ',length=1)
                    R=8.314
                    RT=R*T
                    RS=np.exp((deS-deR)*4185.8518/RT)
                    ee=(RS-1)/(RS+1)
                    printer.res(f'立体选择性ee值:\n{ee*100:.2f}%')
                    
    def writePage(self):
        while True:
            opts=[
                ['1','gif 文件'],
                ['2','xyz 文件'],
                ['3','cub 文件'],
                ['4','si  文件'],
            ]
            opt=self.input.Option('导出文件',opts)
            match opt:
                case None:return
                case '1':
                    from pywfn.writer import gjfWriter
                    for path in self.paths:
                        mol=Mol(get_reader(path))
                        gjfWriter(mol).save()
                case '2':
                    from pywfn.writer import xyzWriter
                    for path in self.paths:
                        mol=Mol(get_reader(path))
                        xyzWriter(mol).save()
                case '3':
                    from pywfn.writer import cubWriter
                    obts=self.input.Number(tip='请输入轨道序号(,)[]: ',type_='int')
                    obts=[obt-1 for obt in obts]
                    step=self.input.Number(tip=f'请输入渲染间隔()[{config.RENDER_CLOUD_STEP}]: ',type_='float',length=1)
                    atoms=self.input.Number(tip='请输入渲染原子(,)[全选]: ',type_='int')
                    for path in self.paths:
                        mol=Mol(reader=get_reader(path))
                        writer=cubWriter(mol)
                        if step is not None:writer.step=step
                        if atoms is not None:writer.atoms=atoms
                        writer.obts=obts
                        writer.save()
                case '4':
                    from pywfn.writer import siWriter
                    printer.info('1.坐标 2.能量 3.频率')
                    selects=self.input.Number(tip='选择导出信息[全选]: ')
                    same=self.input.Bool(tip='保存同一文件?([y]/n): ')
                    for path in printer.track(self.paths,tip='正在导出'):
                        reader=get_reader(path)
                        writer=siWriter(Mol(reader=reader))
                        if selects is not None:writer.selects=selects
                        writer.sameFile=same
                        writer.save()

    def showFiles(self):
        printer.info('当前的读取的文件有：')
        for i,path in enumerate(self.paths):
            print(f'{i},{path.stem}')
        printer.print('')

    def setingPage(self):
        """设置页面"""
        from pywfn.data import temps
        opts=[
            ['1','gif 模板'],
            ['2','si  模板'],
        ]
        while True:
            opt=self.input.Option('设置页面',opts)
            match opt:
                case None:return
                case '1':
                    path,text=self.input.Text()
                    temps.gjf=text
                    config.set_config('TEMPLATE_PATH_GJF',path)
                case '2':
                    path,text=self.input.Text()
                    temps.si=text
                    config.set_config('TEMPLATE_PATH_SI',path)
            printer.res('设置成功')
                
class Input:
    def __init__(self) -> None:
        pass

    def input(self,tip:str):
        res=input(tip)
        if res=='q':sys.exit()
        return res

    def Number(self,tip:str,type_:str='float',length:int=None)->int|float:
        """输入整数"""
        while True:
            opt=input(tip)
            if opt== '':
                return None
            else:
                if not re.match(r'^-?\d+[\.\d+]*([,，]-?\d+[\.\d+]*)*$',opt):
                    printer.warn('格式不正确！！')
                    continue
                nums=re.split(r'[,，]',opt)
                nums=[int(num) if type_=='int' else float(num) for num in nums]
                if length is None:
                    return nums
                else:
                    if len(nums)!=length:
                        continue
                    if length==1:
                        return nums[0]
                    else:
                        return nums
    
    def Bool(self,tip:str,default:bool=True):
        opt=input(tip)
        match opt:
            case '':
                return default
            case 'y':
                return True
            case 'n':
                return False
            case _:
                printer.warn('请输入正确选项!')
                return self.Bool(tip,default)
    
    def Option(self,title,opts,canNone=True):
        # maxLen=max([len(text) for idx,text in opts])
        printer.options(title,opts)
        
        
        idxs=[idx for idx,text in opts]
        while True:
            opt=self.input('请输入对应选项: ')
            if opt=='' and canNone:
                return None
            elif opt in idxs:
                return opt
            else:
                printer.warn('输入选项不存在!')
    
    def Files(self):
        types=['.log','.out','.fch','.gjf'] #支持的文件类型
        while True:
            path=self.input('输入文件[夹]名: ')
            if path=='':continue
            
            path=Path(path)
            if not path.exists():
                printer.wrong('路径不存在')
                continue #如果路径存在
            if path.is_file(): # 如果是文件
                if path.suffix in types:
                    return [path]
                else:
                    printer.warn('不支持的文件类型')
            elif path.is_dir(): # 如果是文件夹
                paths=[]
                for each in path.iterdir(): # 对文件夹中的每个文件进行循环
                    if each.suffix in types: # 如果文件类型是支持的文件类型
                        paths.append(each)
                printer.info(f'共{len(paths)}个文件')
                return paths
    
    def Text(self)->str:
        path=input('请输入文件: ')
        if not Path(path).exists():return None
        return path,Path(path).read_text()