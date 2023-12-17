"""
将多个link的log文件拆分出来
"""
from pathlib import Path
from pywfn.utils import printer
class Tool:
    def __init__(self,path:str) -> None:
        self.path=Path(path)
        self.texts=self.path.read_text()
        self.lines=self.texts.splitlines(keepends=False)
        self.mark=' Initial command:'
        self.nums=[]
        for i,line in enumerate(self.lines):
            if line!=self.mark:continue
            self.nums.append(i)
        self.nums.append(len(self.lines))
    
    def split(self):
        for i in range(len(self.nums)-1):
            l1,l2=self.nums[i:i+2]
            text='\n'.join(self.lines[l1:l2])
            self.write(text,i)
        printer.res('导出完成 >_<')

    def write(self,text,idx):
        name=self.path.stem
        (self.path.parent/f'{name}_{idx+1}.log').write_text(text)
        