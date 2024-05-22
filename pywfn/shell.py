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
        self.mols: dict[str,Mol] = {} # 路径与分子的对应
        printer.ifShell = True
        printer.info(data.start)
        self.input = Inputer(self)

    def homePage(self):
        opts = {
            "0": "导入文件",
            "1": "键の属性",
            "2": "原子属性",
            "3": "实用工具",
            "4": "导出文件",
        }
        while True:
            printer.options('首页', opts)
            opt = input("输入功能选项: ")
            match opt:
                case "0":
                    path=input('请输入文件路径: ')
                    self.input.Files(path)
                case "1":
                    from pywfn import bondprop
                    bondprop.onShell(self)
                case "2":
                    from pywfn import atomprop
                    atomprop.onShell(self)
                case "3":
                    from pywfn import tools
                    tools.onShell(self)
                case "4":
                    from pywfn import writer
                    writer.onShell(self)
                case 'q':
                    break
                case _:
                    printer.warn("命令不存在")


class Inputer:
    def __init__(self, shell: Shell) -> None:
        self.shell = shell

    def input(self, tip: str):
        res = input(tip)
        if res == "q":
            sys.exit()
        return res

    def Number(self, tip: str, dtype: str = "float", length: int = None) -> int | float:
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
                nums = [int(num) if dtype == "int" else float(num) for num in nums]
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

    def Files(self, path:str):
        types = [".log", ".out", ".fch", ".gjf"]  # 支持的文件类型
        pathObj = Path(path)
        if pathObj.is_file():  # 如果是文件
            if pathObj.suffix in types:
                printer.info(f"共1个文件")
                self.shell.paths += [path]
            else:
                printer.warn("不支持的文件类型")
        elif pathObj.is_dir():  # 如果是文件夹
            paths = []
            for each in pathObj.iterdir():  # 对文件夹中的每个文件进行循环
                if each.suffix in types:  # 如果文件类型是支持的文件类型
                    paths.append(each)
            printer.info(f"输入{len(paths)}个文件")
            self.shell.paths += paths

    def Text(self) -> str:
        path = input("请输入文件: ")
        if not Path(path).exists():
            return None
        return path, Path(path).read_text()

    def Moles(self,num:int=None) -> list[Mol]:
        """
        获取当前文件的分子，列举出当前读取的文件让用户选择，将用户选择的文件对应为分子
        """
        paths = self.shell.paths
        printer.info(f"共{len(paths)}个文件:")
        printer.bar()
        for p,path in enumerate(paths):
            print(f'{p:>2} {path}')
        if num is None: # 要满足指定数量
            idxs = self.Number("输入分子编号：", dtype="int")
        else:
            while True:
                idxs = self.Number(f"输入{num}个分子：", dtype="int")
                if len(idxs) == num:break

        mols = []
        for idx in idxs:
            path = paths[idx]
            if path not in self.shell.mols.keys():
                self.shell.mols[path] = Mol(get_reader(path))
            mols.append(self.shell.mols[path])
        return mols