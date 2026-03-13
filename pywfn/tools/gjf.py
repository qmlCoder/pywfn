"""
每个函数都是对gjf文件的修改或生成新的gjf，并返回新的gjf的内容
"""


from pathlib import Path
import numpy as np

from pywfn.base.mole import Mole
from pywfn.reader.gjf import GjfReader
from pywfn.writer.gjf import GjfWriter
from pywfn.editor import Editor
from pywfn.datas import temps
from pywfn.datas.elements import elements

class Tool:
    def __init__(self) -> None:
        pass

    def join(self,paths:list[str])->str:
        """将多个gjf文件拼接到一起"""
        texts=[]
        for path in paths:
            if Path(path).suffix!='.gjf':continue
            text=Path(path).read_text()
            texts.append(text)
        text:str='--Link1--\n\n'.join(texts)
        return text
    