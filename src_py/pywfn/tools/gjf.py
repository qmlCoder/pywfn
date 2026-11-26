"""
每个函数都是对gjf文件的修改或生成新的gjf，并返回新的gjf的内容
"""

from pathlib import Path


class Tool:
    def __init__(self) -> None:
        pass

    def join(self, paths: list[str]) -> str:
        """将多个gjf文件拼接到一起

        Args:
            paths: gjf文件路径列表

        Returns:
            str: 拼接后的gjf文件内容
        """
        texts: list[str] = []
        for path in paths:
            if Path(path).suffix != ".gjf":
                continue
            text = Path(path).read_text()
            texts.append(text)
        result: str = "--Link1--\n\n".join(texts)
        return result
