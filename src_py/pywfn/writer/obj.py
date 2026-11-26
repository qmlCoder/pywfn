"""
obj文件生成器
"""

from pywfn import _core


class ObjWriter:
    def __init__(self, verts=None, faces=None):
        verts = verts if verts is not None else []
        faces = faces if faces is not None else []
        self.core = _core.writer.ObjWriter(verts, faces)  # type: ignore

    @property
    def verts(self):
        return self.core.verts

    @verts.setter
    def verts(self, verts):
        self.core.verts = verts

    @property
    def faces(self):
        return self.core.faces

    @faces.setter
    def faces(self, faces):
        self.core.faces = faces

    def save(self, path: str) -> None:
        """保存obj文件

        Args:
            path: 保存路径
        """
        self.core.save(path)
