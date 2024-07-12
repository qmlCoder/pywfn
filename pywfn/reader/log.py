"""
此脚本用来提取高斯输出文件的信息
高斯的输出文件包含迭代信息(结构优化、扫描等)
但是封装成分子对象之后就只有一个信息了
所以迭代信息只能在reader对象中存在,且不是默认属性
"""

import re
import numpy as np
from pathlib import Path
from functools import lru_cache
from collections import namedtuple
import threading
import json

from pywfn.data.basis import Basis
from pywfn.maths.gto import Gto
from pywfn.utils import printer
from pywfn.data.elements import elements
from pywfn import reader

from typing import Callable
import linecache
import textwrap

class Title:
    def __init__(self,mark:str,jtype:int=0,multi:bool=False) -> None:
        self.line:int=-1
        self.mark:str=mark
        self.jtype=jtype # 0：正则表达式匹配，1：包含式匹配，2：全等匹配
        self.multi=multi # 是否持续查找（多个中查找最后一个）
    
    def set_line(self,line:int):
        self.line=line
    
    def judge(self,line:str):
        """判断所给行是否满足条件"""
        if self.jtype==0:
            return re.search(self.mark,line) is not None
        elif self.jtype==1:
            return self.mark in line
        elif self.jtype==2:
            return line==self.mark

    
    def __repr__(self) -> str:
        return f'{self.line}: {self.mark}'

class LogReader(reader.Reader):
    def __init__(self, path:str,clear:bool=False) -> None:
        super().__init__(path,clear)
        assert type(path)==str,'路径应该为字符串格式'
        assert path[-4:] in ['.log','.out'],'文件类型不匹配，应为.log文件或.out文件'
        self.index=0 #从第一行向下搜索，搜索到需要的就停止
        self.titles={ #记录每一个title所在的行数
            'coords':Title(r'(Input orientation|Standard orientation)',0,True),
            'basis':Title('Standard basis:',1),
            'coefs':Title(r'  Molecular Orbital Coefficients',0),
            'acoefs':Title(r'Alpha Molecular Orbital Coefficients',1),
            'bcoefs':Title(r'Beta Molecular Orbital Coefficients',1),
            'overlap':Title(' *** Overlap *** \n',2),
            'kinetic':Title(' *** Kinetic Energy *** \n',2),
            'potential':Title(' ***** Potential Energy ***** \n',2),
            'density':Title(r'Density Matrix:',1),
            'basisData':Title(r'Overlap normalization',1),
            'engs':Title(r'Zero-point correction',1),
            'keyWards':Title(r'# .+',0),
        }
        
        cpath=Path(f'{self.dfold}/log.json') # config path
        if not cpath.exists():cpath.write_text('{}')
        self.conf:dict=json.loads(cpath.read_text())
        self.cpath=cpath
        self.get_confTitle()
        self.search_title()
    
    def get_confTitle(self):
        """从配置文件中回去标题所在的行"""
        keys=self.conf.keys()
        for tip in self.titles.keys():
            if not f'title_{tip}' in keys:continue
            line=self.conf[f'title_{tip}']
            self.titles[tip].set_line(line)
    def save_config(self):
        self.cpath.write_text(json.dumps(self.conf,indent=4))
    
    def search_title(self):
        """
        在文件中搜索需要的行号
        """
        from concurrent.futures import ThreadPoolExecutor
        from concurrent.futures._base import Future
        # 所有标题所在的行
        def sear_group(start:int): #每一个搜索的线程
            keys=list(self.titles.keys()) # 总的key
            keys=[key for key in keys if self.titles[key].line==-1] # 需要搜索的key
            if len(keys)==0:return
            for j in range(start,start+bsize):
                line=self.getline(j)
                if line=='':break
                for key in keys:
                    title:Title=self.titles[key]
                    if title.line!=-1 and title.multi==False:continue
                    if not title.judge(line):continue
                    self.titles[key].line=j
                    self.conf[f'title_{key}']=j
                    # printer.log(f'搜索到标题,{key},j')
                    self.save_config()
                        # print(j,line)
        nWork=5
        bsize=10_000
        with ThreadPoolExecutor(max_workers=nWork) as executor:
            futures:list[Future] = []
            for i in range(0,self.lineNum,bsize):
                futures.append(executor.submit(sear_group, i))
            for future in futures:
                future.result()
    
    @property
    def normalEnd(self)->bool:
        lastLine=self.getline(self.lineNum)
        return 'Normal termination of Gaussian' in lastLine
    
    @lru_cache
    def get_coords(self)->np.ndarray:
        """原子坐标[n,3]"""
        symbols,coords=self.read_coords()
        assert isinstance(coords,np.ndarray),'coord必须是np.ndarray类型'
        return coords

    @lru_cache
    def get_symbols(self)->list[str]:
        """原子符号[n]"""
        symbols,coords=self.read_coords()
        return symbols

    @lru_cache
    def get_CM(self) -> np.ndarray:
        CM=self.read_CMs()[-1]
        return CM
    
    @lru_cache
    def get_obtAtms(self) -> list[int]:
        ObtAtms=self.read_CMs()[0]
        return ObtAtms
    
    @lru_cache
    def get_obtShls(self) -> list[int]:
        ObtShls=self.read_CMs()[1]
        return ObtShls

    @lru_cache
    def get_obtSyms(self)-> list[list[int]]:
        obtSyms=self.read_CMs()[2]
        return obtSyms

    @lru_cache
    def get_obtEngs(self) -> list[float]:
        ObtEngs=self.read_CMs()[3]
        return ObtEngs
    
    def get_obtOccs(self) -> list[bool]:
        ObtOccs=self.read_CMs()[4]
        return ObtOccs

    @lru_cache
    def get_SM(self)->np.ndarray:
        return self.read_SM()
    
    def get_charge(self) -> int|None:
        read=self.read_multiy()
        if read is None:
            return None
        charge,spin=read
        return charge
    
    def get_spin(self)->int|None:
        read=self.read_multiy()
        if read is None:
            return None
        charge,spin=read
        return spin
    
    def get_energy(self)->float:
        return self.read_energy()
    
    from pywfn.data import Basis
    def get_basis(self)->Basis:
        name=self.read_basisName()
        data=self.read_basisData()
        basis=Basis(name)
        basis.setData(data)
        return basis

    def read_keyWrds(self):
        """读取关键字"""
        titleNum=self.titles['keyWards'].line
        keyWards=''
        for i in range(titleNum,titleNum+3):
            line=self.getline(i,False).strip()
            if line=='-'*70:break
            keyWards+=line
        return keyWards

    @lru_cache
    def read_coords(self)->tuple[list[str],np.ndarray]|None:
        '''读取原子坐标'''
        
        titleNum=self.titles['coords'].line
        assert titleNum!=-1,"没有坐标对应的行"
        if titleNum is None:
            printer.warn('没有读取到原子坐标')
            return None
        s1=r' +\d+ +(\d+) +\d +(-?\d+.\d{6}) +(-?\d+.\d{6}) +(-?\d+.\d{6})'
        coords=[]
        symbols=[]
        for i in range(titleNum+5,self.lineNum):
            line=self.getline(i)
            find=re.search(s1, line)
            if find is not None:
                res=list(find.groups())
                atomID=int(res[0])
                symbol=elements[atomID].symbol
                coord=[float(each) for each in res[1:]]
                symbols.append(symbol)
                coords.append(coord)
            else:
                return symbols,np.array(coords)*1.889

    @lru_cache
    def read_multiy(self)->tuple[int,int]|None:
        """读取电荷和自选多重度"""
        res=re.findall(rf'Charge = +(-?\d) Multiplicity = (\d)',self.text)
        if res is None:
            printer.warn(f'{self.path} 没有读到电荷和自选多重度!')
        else:
            charge,multiy=res[0]
            return int(charge),int(multiy)

    @lru_cache 
    def read_basisName(self)->str:
        """读取基组名"""
        titleNum=self.titles['basis'].line
        if titleNum is None:
            printer.warn('未读取到基组名')
            name='unll'
        else:
            find=re.match(rf' Standard basis: (\S+) ',self.getline(titleNum))
            if find is None:
                name='null'
            else:
                name=find.group(1)
        return name
    
    @lru_cache
    def read_basisData(self):
        """
        读取基组数据
        原子类型
            电子层
        """
        from pywfn.data.basis import BasisData
        assert 'gfinput' in self.read_keyWrds(),'关键词应该包含：gfinput'
        titleNum=self.titles['basisData'].line
        symbols=self.get_symbols()
        if titleNum is None:return
        basisDatas:list[BasisData]=[]
        ifRead=True
        angDict={'S':0,'P':1,'D':2} #角动量对应的字典
        s1='^ +(\d+) +\d+$'
        s2=r' ([SPD]+) +(\d+) \d.\d{2} +\d.\d{12}'
        s3=r'^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})'
        s4=' ****\n'
        atomics=[] # 已经获取过的元素
        atomic=None
        shell=0
        angs=None
        exp=None
        coe=None
        for i in range(titleNum+1,self.lineNum):
            line=self.getline(i)
            if re.search(s1,line) is not None:
                idx=re.search(s1,line).groups()[0] #第几个原子
                idx=int(idx)-1
                symbol=symbols[idx]
                atomic=elements[symbol].idx

                if atomic not in atomics:
                    atomics.append(atomic) #shells
                    ifRead=True # 该元素是否已经读过
                else:
                    ifRead=False
            elif re.search(s2,line) is not None:
                if not ifRead:continue
                shellName,lineNum=re.search(s2,line).groups()
                angs=[angDict[s] for s in shellName] #角动量
                shell+=1
            elif re.search(s3,line) is not None:
                if not ifRead:continue
                numsStr=re.search(s3,line).groups()[0]
                nums:list[str]=re.findall(r'-?\d.\d{10}D[+-]\d{2}',numsStr)
                nums=[float(num.replace('D','E')) for num in nums]
                if len(angs)==1:
                    exp,coe=nums
                    data=BasisData(atomic,shell,angs[0],exp,coe)
                    basisDatas.append(data)
                if len(angs)==2:
                    exp,coe1,coe2=nums
                    data1=BasisData(atomic,shell,angs[0],exp,coe1)
                    data2=BasisData(atomic,shell,angs[1],exp,coe2)
                    basisDatas.append(data1)
                    basisDatas.append(data2)
            elif line==s4 is not None:
                shell=0
                continue
            else:
                break
                # 排序一下
        basisDatas.sort(key=lambda b:(b.atmic,b.shl,b.ang))
        return basisDatas
        
    def read_summery(self):
        '''读取总结信息'''
        summerys=re.findall(r' \d[\\||]\d[\S\s]*?@',self.text)
        summery=summerys[-1]
        summery=''.join([each[1:] for each in summery.split('\n')])
        summeryList=summery.replace('\n','').replace('\\\\','||').split('||')
        Infos,KeyWords,JobTitle,Coords=summeryList[:4]
        _,_,_,jobType,method,basisSet,molFormula,user,date,_=Infos.replace('\\','|').split('|')
        summery={
            'basisSet':basisSet,
            'Coords':np.array([[float(num.replace(' ','')) for num in each.split(',')[1:]] for each in Coords.replace('\\','|').split('|')[1:]])
        }
        return summery
           
    # 分子轨道的文本数据有5种情况
    #情况1                           1         2         3         4         5
    #情况2                           O         O         O         O         O
    #情况3     Eigenvalues --   -11.17917 -11.17907 -11.17829 -11.17818 -11.17794
    #情况4   1 1   C  1S         -0.00901  -0.01132   0.00047  -0.01645  -0.02767
    #情况5   2        2S         -0.00131  -0.00175  -0.00041  -0.00184  -0.00173
    # @lru_cache
    # def read_CM_(self, title:str):  # 提取所有原子的轨道 自己写的代码自己看不懂真实一件可悲的事情,此函数逻辑复杂，要好好整明白
    #     s1=r'^(( +\d+){1,5}) *$'
    #     s2=r'^(( *(\(\w+\)--)?[OV]){1,5})$'
    #     s3=r'^ +Eigenvalues --(( +-?\d+.\d+){1,5})'
    #     s4=r'^ +\d+ +(\d+) +([A-Za-z]+) +(\d[A-Z]+)(( *-?\d+.\d+){1,5})$'
    #     # s5='^ +\d+ +(\d+[A-Za-z ]+)(( *-?\d+.\d+){1,5})$'
    #     s5=r'^ +\d+ +(\d+[A-Za-z]+ ?\+?-?\d?)(( *-?\d+.\d+){1,5})$'

    #     titleNum=self.titles[title].line
    #     if titleNum is None:return None

    #     OrbitalAtom=[] # 系数矩阵每一行对应的原子
    #     OrbitalLayer=[] # 系数矩阵每一行对应的壳层符号
    #     OrbitalEngs=[] # 轨道本征值
    #     OrbitalType=[] # 轨道类型
    #     dataDict:dict[int:list[list[float]]]={}
    #     firstShow=True
    #     if titleNum is None:
    #         return
    #     for i in range(titleNum+1,self.lineNum):
    #         line=self.getline(i)
    #         if re.search(s1, line) is not None: #情况1
    #             pass
    #         elif re.search(s2, line) is not None: # 情况2，获得column
    #             OrbitalType += re.split(r' +', line.replace('\n',''))[1:] # 获取占据轨道还是非占据轨道
                
    #         elif re.search(s3, line) is not None: # 情况3，每个轨道本征值     
    #             line_data,_=re.search(s3,line).groups()
    #             line_data=re.findall('-?\d+.\d+', line_data)
    #             line_data=[float(each) for each in line_data]
    #             OrbitalEngs+=line_data
    #         elif re.search(s4,line) is not None: # 原子出现
    #             atomIDX,atomType,layer,line_data,_=re.search(s4,line).groups()
    #             atomIDX = int(atomIDX)
    #             nums=[float(each) for each in re.findall('-?\d+.\d+',line_data)]
    #             if atomIDX not in dataDict.keys():
    #                 dataDict[atomIDX]=[]
    #             else:
    #                 firstShow=False
    #             dataDict[atomIDX].append([])
    #             dataDict[atomIDX][-1].append(nums)
    #             if firstShow:OrbitalLayer.append(layer)
    #             if firstShow:OrbitalAtom.append(atomIDX)
    #             # self.OCdict[atomIDX].set(layer,nums)
    #         elif re.search(s5,line) is not None:
    #             layer,line_data,_=re.search(s5,line).groups()
    #             line_data=re.findall('-?\d+.\d+', line_data)
    #             nums=[float(each) for each in line_data]
    #             dataDict[atomIDX][-1].append(nums)
    #             if firstShow:OrbitalLayer.append(layer)
    #             if firstShow:OrbitalAtom.append(atomIDX)
    #             # self.OCdict[atomIDX].set(layer,nums)
    #         else: # 若不满足以上任意一种情况，说明已经查找完毕，则对收集到的数据进行处理
    #             printer.console.log(f'读取完成,i={i},line={line}')
    #             for atomic,matrics in dataDict.items():
    #                 for i,matrix in enumerate(matrics):
    #                     dataDict[atomic][i]=np.array(matrix)
    #             for atomic,matrics in dataDict.items():
    #                 dataDict[atomic]=np.concatenate(matrics,axis=1)
    #             CM=np.concatenate([m for m in dataDict.values()],axis=0)
    #             break
    #     printer.console.log(f'读取完成CM,shape={CM.shape}')
    #     assert CM.shape[0]==CM.shape[1],"CM需要为正方形矩阵"
    #     return OrbitalAtom,OrbitalLayer,OrbitalEngs,OrbitalType,CM
    
    
    @lru_cache
    def read_CM(self, title:str)->tuple[list[int],list[int],list[str],list[float],list[bool],np.ndarray]:  # 提取所有原子的轨道 自己写的代码自己看不懂真实一件可悲的事情,此函数逻辑复杂，要好好整明白
        
        keyWards=self.read_keyWrds()
        assert 'pop=full' in keyWards,'关键词应包含：pop=full'
        find=re.search(r'NBasis *= *(\d+)',self.text)
        assert find is not None,"没有找到基函数数量!!"
        NBasis=int(find.groups()[0])
        NBlock=NBasis//5+(0 if NBasis%5==0 else 1)
        titleNum=self.titles[title].line
        
        assert titleNum!=-1,"没有轨道系数!!"
        blockLen=NBasis+3 #一块数据行数，在log文件的输出中，轨道系数是按照每一块五列来分块的
        # print(titleNum,NBlock,blockLen)
        ObtAtms=[] # 原子轨道对应的原子
        ObtShls=[] # 占据原子的第多少层轨道
        ObtSyms=[] # 原子轨道对应的角动量
        ObtEngs=[] # 分子轨道能量
        ObtOccs=[] # 占据/非占据
        
        CM=np.zeros(shape=(NBasis,NBasis))
        for i,l in enumerate(range(titleNum+1,titleNum+NBlock*blockLen+1,blockLen)):
            ltex=self.getline(l+1)
            occs=self.getline(l+1)[21:71] # 占据类型的起止位置
            occs=[occs[i:i+10].strip() for i in range(0,50,10)]
            occs=[e for e in occs if e!='']
            occs=[True if e[-1]=='O' else False for e in occs]
            ObtOccs+=occs
            # engs =re.split(' +',self.getline(l+2)[21:].strip())
            engs =textwrap.wrap(self.getline(l+2)[21:].strip(),10)
            line=self.getline(l+2)
            engs=[float(e) for e in engs]
            ObtEngs+=engs
            
            coefs=[line[21:] for line in self.getlines(l+3,l+3+NBasis)]
            for j,line in enumerate(coefs):
                coef=re.findall(r'-?\d+.\d+',line)
                coef=[float(e) for e in coef]
                CM[j,i*5:i*5+len(coef)]=coef

            if i==0: #只在第一块记录行信息
                atomID=''
                for l2 in range(l+3,l+blockLen):
                    line=self.getline(l2)[:16] # 行信息(角动量，对应原子)的起止位置
                    match=line[12:].strip()
                    shl,sym=match[0],match[1:] #这里壳层只能是整数哦
                    ObtShls.append(int(shl))
                    ObtSyms.append(sym)
                    obtAtom=line[5:9].strip()
                    if obtAtom!='':
                        atomID=int(obtAtom) # 更改当前行的原子
                    ObtAtms.append(atomID)
        return ObtAtms,ObtShls,ObtSyms,ObtEngs,ObtOccs,CM

    @lru_cache
    def read_CMs(self)->tuple[list,list,list,list,list,np.ndarray]:
        """
        获取轨道系数及相关信息
        atms,shls,angs,engs,occs,CM
        """
        atms=self.load_fdata('atms.npy')
        shls=self.load_fdata('shls.npy')
        angs=self.load_fdata('angs.npy')
        engs=self.load_fdata('engs.npy')
        occs=self.load_fdata('occs.npy')
        CM=self.load_fdata('CM.npy')
        fdatas=[atms,shls,angs,engs,occs,CM]
        rdatas=[fdata is not None for fdata in fdatas]
        if all(rdatas): # 所有读取的数据都不为空，有可能某些数据被用户删除
            lists:list[np.ndarray]=[atms,shls,angs,engs,occs]
            atms,shls,angs,engs,occs=[e.tolist() for e in lists]
            return atms,shls,angs,engs,occs,CM
        elif self.titles['coefs'].line!=-1:
            atms,shls,angs,engs,occs,CM=self.read_CM('coefs')
        elif self.titles['acoefs'].line!=-1:
            atmsA,shlsA,angsA,engsA,occsA,CMA=self.read_CM('acoefs')
            atmsB,shlsB,angsB,engsB,occsB,CMB=self.read_CM('bcoefs')
            atms=atmsA
            angs=angsA
            shls=shlsA
            engs=engsA+engsB
            occs=occsA+occsB
            CM=np.concatenate([CMA,CMB],axis=1)
        else:
            raise ValueError('没有找到轨道系数')
        self.save_fdata('atms',atms)
        self.save_fdata('shls',shls)
        self.save_fdata('angs',angs)
        self.save_fdata('engs',engs)
        self.save_fdata('occs',occs)
        self.save_fdata('CM',CM)
        return atms,shls,angs,engs,occs,CM

    @lru_cache
    def read_SM(self):
        """读取重叠矩阵"""
        SM=self.load_fdata('SM.npy')
        if SM is None:
            SM=self.read_Mat('overlap')
            self.save_fdata('SM.npy',SM)
        return SM

    def read_Mat(self,title:str)->np.ndarray:
        """读取矩阵"""
        s1=r'^( +\d+){1,5} *$'
        s2=r' +\d+( +-?\d.\d{6}D[+-]\d{2}){1,5} *'
        lineDatas:dict[str,str]={}
        titleNum=self.titles[title].line
        for i in range(titleNum+1,self.lineNum):
            line=self.getline(i)
            if re.match(s1, line) is not None:
                pass
            elif re.match(s2, line) is not None:
                idx=line[:8]
                if idx not in lineDatas.keys():
                    lineDatas[idx]=''
                nums=line[8:].replace('D','e')
                lineDatas[idx]+=f'{nums} '
            else:
                break
        size=len(lineDatas.keys())
        SM_=np.zeros((size,size))
        for idx,each in lineDatas.items():
            each:str=each
            nums=re.split(' +',each.strip())
            nums=[float(num) for num in nums]
            SM_[int(idx)-1,:len(nums)]=nums
        return np.tril(SM_)+np.tril(SM_,-1).T

    @lru_cache
    def read_charge(self)->float:
        res=re.search(r'Sum of Mulliken charges = +(-?\d.\d{5})',self.text).group(1)
        return float(res)
    
    def read_energy(self)->float:
        engs=re.findall(r'SCF Done: +E\(.*\) += +(-?\d+.\d+)',self.text)
        if engs:return float(engs[-1])
        printer.warn('未读取到能量，使用0代替')
        return 0.0

    @lru_cache
    def read_energys(self)->list[float]:
        """返回表格中对应的这8个能量，只返回数值"""
        engList = [
            'Zero-point correction', 
            'Thermal correction to Energy', 
            'Thermal correction to Enthalpy',
            'Thermal correction to Gibbs Free Energy', 
            'Sum of electronic and zero-point Energies',
            'Sum of electronic and thermal Energies', 
            'Sum of electronic and thermal Enthalpies',
            'Sum of electronic and thermal Free Energies'
            ]
        searhNum=0
        engDict={e:0.0 for e in engList}
        lineNum=self.titles['engs'].line
        if lineNum is None:printer.log("{self.path} 未读取到能量")
        for i in range(lineNum,self.lineNum):
            line=self.getline(i)
            for each in engList:
                res=re.search(rf'{each}=\s+(-?\d+\.\d+)',line)
                if res is not None:
                    re.findall(r'-?\d+\.\d+',line)
                    engDict[each]=float(res.groups()[0])
                    searhNum+=1
            if searhNum==len(engList):break
        engNums=[float(e) for e in engDict.values()]
        return engNums
    
    def read_freqs(self):
        """读取频率"""
        freqs=[]
        find = re.findall(r'^\s+Frequencies\s--\s+(-?\d+\.\d+.+)$', self.text, flags=re.M)
        if find is None:return freqs
        for line in find:
            finds=re.findall('-?\d+.\d+',line)
            freqs+=[float(e) for e in finds]
        return freqs
        
