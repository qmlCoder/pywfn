from pywfn import core


class Tool:
    def __init__(self) -> None:
        self.core=core.tools.log.Tool() # type: ignore

    def get_SI(self,paths:list[str]) -> str:
        return self.core.get_SI(paths)
    
    def split_opt(self,path:str):
        return self.core.split_opt(path)
    
    # def split_scan(self,path:str) -> list[Mole]:
    #     start=0
    #     read=False
    #     geomes:list[Geome]=[]
    #     engs:list[float]=[]
    #     with open(path,'r') as file:
    #         for l,line in enumerate(file):
    #             line=line.rstrip()
    #             # print(line)
    #             if line.strip()=='Input orientation:':
    #                 read=True
    #                 start=l+5
    #                 geomes.append(Geome())
    #             if l>=start and line.strip()=='-'*69:
    #                 read=False
    #             if read and l>=start:
    #                 match=re.match(r" +\d+ +(\d+) +\d+ +(-?\d+.\d+) +(-?\d+.\d+) +(-?\d+.\d+)",line)
    #                 assert match is not None,"匹配结果不正确"
    #                 atomic,x,y,z=match.groups()
    #                 sym=elements[int(atomic)].sym
    #                 x=float(x)
    #                 y=float(y)
    #                 z=float(z)
    #                 # print(sym,x,y,z)
    #                 geomes[-1].addAtom(sym,np.array([x,y,z]))
    #             if line[:23]==' SCF Done:  E(RB3LYP) =':
    #                 # print(line)
    #                 eng=re.search(r'-?\d+.\d+',line)
    #                 assert eng is not None,"未匹配到能量"
    #                 eng=eng.group()
    #                 engs.append(float(eng))
    #             if line==' Optimization completed.':
    #                 break
    #     moles:list[Mole]=[]
    #     for g,geome in enumerate(geomes):
    #         reader=NonReader()
    #         reader.geome=geome
    #         mol=Mole(reader)
    #         mol._props['energy']=engs[g]
    #         moles.append(mol)
    #     return moles

    # def split_irc(self):
    #     pass

    # def split_link(self,path:str):
    #     """将link链接的作业分割成单独的log文件"""
    #     def save(path:str,text:str):
    #         with open(path,'w',encoding='utf-8') as f:
    #             f.write(text)
    #         print(f'save {path} complete')      
    #     root=Path(path).parent
    #     stem=Path(path).stem
    #     self.fold=root/f'{stem}_lsp'
    #     self.fold.mkdir(exist_ok=True)
    #     lkIdx=0
    #     ftext=''
    #     with open(path,'r',encoding='utf-8') as file:
    #         for line in file:
    #             ftext+=f'{line}'
    #             if 'Normal termination of Gaussian' in line: # 如果遇到标志行
    #                 save(f'{self.fold}/{lkIdx}.log',ftext)
    #                 lkIdx+=1
    #                 ftext=''
    #     print('导出完成 >_<')
            