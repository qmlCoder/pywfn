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
from pywfn import tools
from pywfn import config
from pywfn.utils import printer
from pywfn import data
from pywfn.base import Mol
import numpy as np


class Shell:
    def __init__(self):
        self.paths: list[Path] = []
        printer.ifShell = True
        printer.info(data.start)
        self.input = Input(self)

    def homePage(self):
        opts = {
            "1": "键属性",
            "2": "原子属性",
            "3": "实用工具",
            "4": "导出文件",
            "5": "程序设置",
        }
        while True:
            opt = self.input.Option("首页", opts, must=True)
            match opt:
                case "1":
                    self.bondProp()
                case "2":
                    self.atomProp()
                case "3":
                    self.toolsPage()
                case "4":
                    self.writePage()
                case "5":
                    self.setingPage()
                case _:
                    printer.warn("命令不存在")

    def atomProp(self):
        """计算原子属性"""
        opts = {
            "1": "电荷分布",
            "2": "电子自旋",
            "3": "π 电子性质",
            "4": "原子自由价",
            "5": "方向电子性质",
            "6": "福井/parr函数",
        }
        chrgs={'1':'mulliken','2':'lowdin'}
        props={'1':'charge','2':'spin'}
        funcs={'1':'fukui','2':'parr'}
        while True:
            opt = self.input.Option("计算原子属性", opts)
            match opt:
                case None:
                    break
                case "1":  # 电荷分布
                    from pywfn.atomprop import atomCharge
                    chrg = self.input.Option(title="输入电荷类型", opts=chrgs, must=True)
                    mol = self.input.Moles()[0]
                    caler = atomCharge.Calculator(mol)
                    if chrg == "1":
                        caler.chrg = "mulliken"
                    if chrg == "2":
                        caler.chrg = "lowdin"
                    caler.logRes()

                case "2":  # 电子自旋
                    chrgs = [["1", "mulliken"], ["2", "lowdin"]]
                    chrg = self.input.Option(title="输入自旋类型", opts=chrgs, must=True)
                    from pywfn.atomprop import atomSpin

                    mol = self.input.Moles()[0]
                    caler = atomSpin.Calculator(mol)
                    if chrg == "1":
                        caler.chrg = "mulliken"
                    if chrg == "2":
                        caler.chrg = "lowdin"
                    caler.logRes()

                case "3":  # π 电子性质
                    from pywfn.atomprop import piProps

                    mol = self.input.Moles()[0]
                    caler = piProps.Calculator(mol)
                    caler.chrg = self.input.Option(title="电荷类型", opts=chrgs, must=True)
                    caler.prop = self.input.Option(title='计算性质', opts=props, must=True)
                    caler.print(caler.resStr())

                case "4":  # 原子自由价
                    from pywfn.atomprop import freeValence

                    mol = self.input.Moles()[0]
                    caler = freeValence.Calculator(mol)
                    while True:
                        idx = self.input.Number(
                            tip="输入原子编号: ", type_="int", length=1
                        )
                        if idx is None:
                            break
                        caler.vect = self.input.Vector(tip="输入法向量", norm=True)
                        caler.logRes()

                case "5":  # 方向电子性质
                    from pywfn.atomprop import dirProps
                    
                    mol = self.input.Moles()[0]
                    caler = dirProps.Calculator(mol)
                    caler.chrg = self.input.Option(title="电荷类型", opts=chrgs, must=True)
                    caler.prop = self.input.Option(title='计算性质', opts=props, must=True)
                    
                    
                    while True:
                        atom = self.input.Number(
                            tip="输入原子编号: ", type_="int", length=1
                        )
                        if atom is None:
                            break
                        vect = self.input.Vector(tip="输入方向向量")
                        if vect is None:
                            break
                        caler.atoms = [atom]
                        caler.vects = [vect]
                        caler.printRes()

                case "6":  # 福井函数、parr函数
                    from pywfn.atomprop import delProps

                    printer.info("分别输入-,0,+对应的分子编号：")
                    mols = self.input.Moles()
                    caler = delProps.Calculator(mols)
                    caler.chrg = self.input.Option(title="电荷类型", opts=chrgs, must=True)
                    caler.prop = self.input.Option(title='计算性质', opts=props, must=True)
                    caler.func = self.input.Option(title='使用函数', opts=funcs, must=True)
                    while True:
                        atom = self.input.Number(
                            tip="输入原子编号: ", type_="int", length=1
                        )
                        if atom is None:
                            break
                        vect = self.input.Vector(tip="输入方向向量：", norm=True)
                        if vect is None:
                            break
                        caler.atoms = [atom]
                        caler.vects = [vect]
                        ev,nv=caler.calculate().flatten() # [1,2]
                        printer.res(f"方向亲电/亲核指标：{ev:.4f},{nv:.4f}")
                        
                case "7":  # 获取各种方向
                    while True:
                        idx = self.input.Number(
                            tip="输入原子编号: ", type_="int", length=1
                        )
                        if idx is None:
                            break
                        atom = mol.atom(idx)

                        normal = atom.get_Normal()
                        obtWay = atom.get_obtWay(len(mol.O_obts) - 1)

                        if normal is not None:
                            printer.vector("原子法向量", normal)
                        else:
                            printer.warn("原子没有法向量")
                        printer.vector("HOMO轨道方向", obtWay)
                case _:
                    return

    def bondProp(self):
        """计算各种键级"""
        opts = {
            "1": "π 键级(OP)",
            "2": "π 键级(SMO)",
            "3": "Mayer 键级",
            "4": "HMO 键级",
        }
        mol=self.input.Moles()[0]
        while True:
            match self.input.Option("计算各种键级", opts):
                case None:
                    break
                case "1":
                    from pywfn.bondprop import piDM

                    caler = piDM.Calculator(mol)
                    while True:
                        bond = self.input.Number(
                            "输入两原子的编号(,)[]: ", type_="int", length=2
                        )
                        if bond is None:
                            break
                        caler.bond = bond
                        caler.print()
                case "2":
                    from pywfn.bondprop import piSM

                    caler = piSM.Calculator(mol)
                    while True:
                        bond = self.input.Number(
                            "输入两原子的编号(,)[]: ", type_="int", length=2
                        )
                        if bond is None:
                            break
                        caler.bond = bond
                        caler.print()
                case "3":
                    from pywfn.bondprop import mayer

                    caler = mayer.Calculator(mol)
                    while True:
                        bond = self.input.Number(
                            "输入两原子的编号(,)[]: ", type_="int", length=2
                        )
                        if bond is None:
                            break
                        caler.bond = bond
                        caler.print()
                case "4":
                    from pywfn.bondprop import hmo

                    caler = hmo.Calculator(mol)
                    while True:
                        bond = self.input.Number(
                            "输入两原子的编号(,)[]: ", type_="int", length=2
                        )
                        if bond is None:
                            break
                        caler.bond = bond
                        caler.print()
                case _:
                    printer.warn("选项不存在")
                    continue

    def toolsPage(self):
        """进入实用工具页面"""
        opts = {
            "1": "分割 扫描 文件",
            "2": "分割 IRC 文件",
            "3": "分割 link 任务",
            "4": "生成 PES 文件",
            "5": "计算 反应活化能",
            "6": "计算 立体选择性ee值",
        }
        while True:
            opt = self.input.Option("使用工具页面", opts)
            match opt:
                case None:
                    break
                case "1":
                    for i, path in enumerate(self.paths):
                        tools.SplitScan(path).split()
                    printer.res("文件分割完成 >_<")
                case "2":
                    for i, path in enumerate(self.paths):
                        tools.SplitIrc(path).split()
                    printer.res("文件分割完成 >_<")
                case "3":
                    for i, path in enumerate(self.paths):
                        tools.SplitLink(path).split()
                case "4":
                    from pywfn.tools.getPES import Tool

                    self.showFiles()
                    tool = Tool(self.paths)
                    tool.route = self.input.Number(tip="路径索引：", type_="int")
                    tool.create()
                case "5":
                    from pywfn.chemProp import reaActEne
                    K = self.input.Number(tip="输入K(催化常数): ", length=1)
                    T = self.input.Number(tip="输入T(反应温度,K): ", length=1)
                    g0 = reaActEne(K,T)
                    g1 = g0 / 4.184
                    printer.res(
                        f"G(吉布斯自由能):\n{g0:>10.4f} kJ/mol\n{g1:>10.4f} Kcal/mol"
                    )
                case "6":
                    from pywfn.chemProp import steSelEE
                    deR, deS = self.input.Number(tip="分别输入ΔG(R),ΔG(S): ", length=2)
                    T = self.input.Number(tip="输入T(反应温度,K): ", length=1)
                    ee = steSelEE(deR,deS,T)
                    printer.res(f"立体选择性ee值:\n{ee*100:.2f}%")

    def writePage(self):
        while True:
            opts = {
                "1": "gif 文件",
                "2": "xyz 文件",
                "3": "cub 文件",
                "4": "si  文件",
            }
            opt = self.input.Option("导出文件", opts)
            match opt:
                case None:
                    return
                case "1":
                    from pywfn.writer import gjfWriter

                    for path in self.paths:
                        mol = Mol(get_reader(path))
                        gjfWriter(mol).save()
                case "2":
                    from pywfn.writer import xyzWriter

                    for path in self.paths:
                        mol = Mol(get_reader(path))
                        xyzWriter(mol).save()
                case "3":
                    from pywfn.writer import cubWriter

                    obts = self.input.Number(tip="请输入轨道序号(,)[]: ", type_="int")
                    obts = [obt - 1 for obt in obts]
                    step = self.input.Number(
                        tip=f"请输入渲染间隔()[{config.RENDER_CLOUD_STEP}]: ",
                        type_="float",
                        length=1,
                    )
                    atoms = self.input.Number(
                        tip="请输入渲染原子(,)[全选]: ", type_="int"
                    )
                    direct = self.input.Vector(
                        tip="请输入投影方向[]", norm=True, must=False
                    )
                    for path in self.paths:
                        mol = Mol(reader=get_reader(path))
                        writer = cubWriter(mol)
                        if step is not None:
                            writer.step = step
                        if atoms is not None:
                            writer.atoms = atoms
                        if direct is not None:
                            writer.direct = direct
                        writer.obts = obts
                        writer.save("")
                case "4":
                    from pywfn.writer import siWriter

                    printer.info("1.坐标 2.能量 3.频率")
                    selects = self.input.Number(tip="选择导出信息[全选]: ")
                    same = self.input.Bool(tip="保存同一文件?([y]/n): ")
                    for path in printer.track(self.paths, tip="正在导出"):
                        reader = get_reader(path)
                        writer = siWriter(Mol(reader=reader))
                        if selects is not None:
                            writer.selects = selects
                        writer.sameFile = same
                        writer.save()

    def showFiles(self):
        printer.info(f"当前读取{len(self.paths)}文件：")
        for i, path in enumerate(self.paths):
            print(f"{i},{path.stem}")
        printer.print("")

    def setingPage(self):
        """设置页面"""
        from pywfn.data import temps

        opts = {
            "1": "gif 模板",
            "2": "si  模板",
        }
        while True:
            opt = self.input.Option("设置页面", opts)
            match opt:
                case None:
                    return
                case "1":
                    path, text = self.input.Text()
                    temps.gjf = text
                    config.set_config("TEMPLATE_PATH_GJF", path)
                case "2":
                    path, text = self.input.Text()
                    temps.si = text
                    config.set_config("TEMPLATE_PATH_SI", path)
            printer.res("设置成功")


class Input:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell

    def input(self, tip: str):
        res = input(tip)
        if res == "q":
            sys.exit()
        return res

    def Number(self, tip: str, type_: str = "float", length: int = None) -> int | float:
        """输入整数"""

        def split1(frg: str):
            if "-" in frg:
                num1, num2 = [int(e) for e in frg.split("-")]
                assert num1 < num2, "数字1要大于数字2"
                return list(f"{e}" for e in range(num1, num2 + 1))
            else:
                return [frg]

        while True:
            opt = input(tip)
            if opt == "":
                return None
            else:
                if not re.match(r"^-?\d+[\.\d+]*([,，-]-?\d+[\.\d+]*)*$", opt):
                    printer.warn("格式不正确！！")
                    continue

                frgs = re.split(r"[,，]", opt)
                nums = []
                for frg in frgs:
                    nums += split1(frg)
                nums = [int(num) if type_ == "int" else float(num) for num in nums]
                if length is None:
                    return nums
                else:
                    if len(nums) != length:
                        continue
                    if length == 1:
                        return nums[0]
                    else:
                        return nums

    def Vector(self, tip: str, norm: bool = False, must=False) -> np.ndarray:
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

    def Bool(self, tip: str, default: bool = True):
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
                return self.Bool(tip, default)

    def Option(self, title: str, opts:dict[str,str], must=False):
        # maxLen=max([len(text) for idx,text in opts])
        printer.options(title, opts)
        idxs = opts.keys()
        while True:
            opt = self.input("请输入对应选项: ")
            if opt == "":
                if must:
                    continue
                else:
                    return None
            elif opt in idxs:
                return opt
            elif Path(opt).exists():  # 如果输入的是一个文件
                self.Files(opt)
            elif opt == "sf":
                self.shell.showFiles()
            elif opt == "cf":
                printer.info("清空文件")
                self.shell.paths = []
                self.shell.showFiles()
            else:
                printer.warn("无效的输入!")

    def Files(self, path):
        types = [".log", ".out", ".fch", ".gjf"]  # 支持的文件类型
        path = Path(path)
        if path.is_file():  # 如果是文件
            if path.suffix in types:
                printer.info(f"共1个文件")
                self.shell.paths += [path]
            else:
                printer.warn("不支持的文件类型")
        elif path.is_dir():  # 如果是文件夹
            paths = []
            for each in path.iterdir():  # 对文件夹中的每个文件进行循环
                if each.suffix in types:  # 如果文件类型是支持的文件类型
                    paths.append(each)
            printer.info(f"输入{len(paths)}个文件")
            self.shell.paths += paths

    def Text(self) -> str:
        path = input("请输入文件: ")
        if not Path(path).exists():
            return None
        return path, Path(path).read_text()

    def Moles(self) -> list[Mol]:
        paths = self.shell.paths
        if len(paths) == 1:
            path = paths[0]
            reader = get_reader(path)
            return [Mol(reader=reader)]
        self.shell.showFiles()
        idxs = self.Number("输入分子编号：", type_="int")

        mols = []
        for idx in idxs:
            path = paths[idx]
            reader = get_reader(path)
            mol = Mol(reader=reader)
            mols.append(mol)
        return mols
