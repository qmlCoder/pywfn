"""
obj文件生成器
"""

from pywfn import _core


class ObjWriter:
    """OBJ 三维模型文件写出器"""

    def __init__(self, verts: list[list[float]] | None = None, faces: list[list[int]] | None = None):
        """obj文件写出器

        Args:
            verts: 顶点坐标列表
            faces: 三角面片列表，每个面片是三个顶点索引
        """
        verts = verts if verts is not None else []
        faces = faces if faces is not None else []
        self.core = _core.writer.ObjWriter(verts, faces)  # type: ignore

    @property
    def verts(self) -> list[list[float]]:
        """顶点坐标列表"""
        return self.core.verts

    @verts.setter
    def verts(self, verts: list[list[float]]) -> None:
        self.core.verts = verts

    @property
    def faces(self) -> list[list[int]]:
        """三角面片列表"""
        return self.core.faces

    @faces.setter
    def faces(self, faces: list[list[int]]) -> None:
        self.core.faces = faces

    def save(self, path: str) -> None:
        """保存obj文件

        Args:
            path: 保存路径
        """
        self.core.save(path)
