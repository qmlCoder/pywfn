"""
定义分子结构编辑器
可以用来生成刚性扫描的结构
"""

from pywfn import _core
from pywfn.base import Mole


class Editor:
    def __init__(self, mole: Mole) -> None:
        self.mole = mole
        self.core = _core.editor.Editor(mole.core)  # type: ignore
