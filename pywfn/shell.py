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
from pywfn import utils
import numpy as np


class Shell:
    def __init__(self):
        self.paths: list[str] = []
        self.mols: dict[str,Mol] = {} # 路径:分子的对应
        printer.ifShell = True
        printer.info(data.start)
        self.input = Inputer(self)

    def homePage(self):
        opts = {
            "0": "导入文件",
            "1": "键の性质",
            "2": "原子性质",
            "3": "实用工具",
            "4": "导出文件",
        }
        while True:
            opts["0"]=f"导入文件({len(self.paths)})"
            printer.options('首页', opts)
            opt = input("输入功能选项: ")
            match opt:
                case "0":
                    self.input.Files()
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
    
    def Float(self,tip:str,count:int=0)->list[float]:
        """输入浮点数"""
        while True:
            numStr=input(tip)
            if numStr=="":continue
            nums=[float(e) for e in numStr.split(',')]
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
            if numStr=="":continue
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
                printer.warn(f"文件数量不够,至少需要{count}个")
                self.Files()
                continue
            idxs = self.Integ(f"输入分子编号: ")
            if len(idxs) == count:break
            if count==0:break
        return [paths[int(idx)] for idx in idxs]

    def Moles(self,count:int=0) -> list[Mol]:
        """
        获取当前文件的分子，列举出当前读取的文件让用户选择，将用户选择的文件对应为分子
        num:分子的数量
        mtype:返回路径还是分子类
        """
        paths=self.Paths(count)

        mols = []
        for path in paths:
            if path not in self.shell.mols.keys():
                self.shell.mols[path] = Mol(get_reader(path))
            mols.append(self.shell.mols[path])
        
        return mols