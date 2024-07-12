"""
将多个gjf拼接到一起
"""
from pywfn.utils import printer

from pathlib import Path

class Tool:
    def __init__(self,paths:list[str]) -> None:
        self.paths=paths

    def build(self):
        texts=[]
        for path in self.paths:
            if Path(path).suffix!='.gjf':continue
            text=Path(path).read_text()
            texts.append(text)
        return texts
    
    def save(self,path:str):
        texts=self.build()
        Path(path).write_text('--Link1--\n\n'.join(texts),encoding='utf-8')
        printer.info(f'文件成功导出至{path} >_<')