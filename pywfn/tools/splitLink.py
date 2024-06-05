"""
将多个link的log文件拆分出来
log->log
文件会很大不要全部读取
根据“ Initial command:”进行分割
每个分隔标志往下到下一个分割标志为一段
"""
from pathlib import Path
from pywfn.utils import printer
class Tool:
    def __init__(self,path:str) -> None:
        self.path=path
        self.mark=' Initial command:\n'
        root=Path(path).parent
        stem=Path(path).stem
        self.fold=root/f'{stem}_lsp'
        if not self.fold.exists():
            self.fold.mkdir()
        self.lkIdx=0
    
    def split(self):
        lkIdx=0
        ftext=''
        with open(self.path,'r',encoding='utf-8') as file:
            for line in file:
                if line==self.mark: # 如果遇到标志行
                    if lkIdx!=0:
                        self.save(f'{self.fold}/{lkIdx}.log',ftext)
                    print(f'spliting {lkIdx} ...')
                    lkIdx+=1
                    ftext=''
                else:
                    ftext+=f'{line}'
            self.save(f'{self.fold}/{lkIdx}.log',ftext)
        printer.res('导出完成 >_<')
    
    def save(self,path,text):
        with open(path,'w',encoding='utf-8') as f:
            f.write(text)
        print(f'save {path} complete')       