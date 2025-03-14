"""
此脚本用来提取高斯输出文件的信息
高斯的输出文件包含迭代信息(结构优化、扫描等)
但是封装成分子对象之后就只有一个信息了
所以迭代信息只能在reader对象中存在,且不是默认属性
输出文件中存储的系数可能是球谐的，但是要转为笛卡尔的才方便使用
"""

import re
import numpy as np
from pathlib import Path
from functools import lru_cache
from collections import namedtuple
import threading
import json

from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.utils import printer
from pywfn.data.elements import elements
from pywfn import reader

from typing import Callable
import linecache
import textwrap
import time
from pywfn.reader.utils import toCart


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
        if self.jtype==0: #正则表达式匹配
            return re.search(self.mark,line) is not None
        elif self.jtype==1: #包含
            return self.mark in line
        elif self.jtype==2: #相等
            return line==self.mark
    
    def set_lnum(self,lnum:int):
        if lnum>=self.line:
            self.line=lnum

    
    def __repr__(self) -> str:
        return f'{self.line}: {self.mark}'

class LogReader(reader.Reader):
    def __init__(self, path:str,cache:bool=False) -> None:
        super().__init__(path,cache)
        assert type(path)==str,'路径应该为字符串格式'
        assert Path(path).suffix in ['.log','.out'],'文件类型不匹配，应为.log文件或.out文件'
        self.index=0 #从第一行向下搜索，搜索到需要的就停止
        self.titles={ #记录每一个title所在的行数
            'coords':Title(r'(Input orientation|Standard orientation)',0,True),
            'basis':Title('Standard basis:',1),
            'coefs':Title(r'  Molecular Orbital Coefficients',0,True),
            'acoefs':Title(r'Alpha Molecular Orbital Coefficients',1,True),
            'bcoefs':Title(r'Beta Molecular Orbital Coefficients',1,True),
            'overlap':Title(' *** Overlap *** \n',2),
            'kinetic':Title(' *** Kinetic Energy *** \n',2),
            'potential':Title(' ***** Potential Energy ***** \n',2),
            'density':Title(r'Density Matrix:',1),
            'basisData':Title(r'Overlap normalization',1),
            'engs':Title(r'Zero-point correction',1),
            'keyWards':Title(r'^ # .+',0),
            'nele':Title(r'alpha electrons',1)
        }
        self.conf={}
        self.cpath=Path(f'{self.dataFold}/log.json') # config path 配置文件，减小第二次使用搜索时间的消耗
        if self.cache:
            if not self.cpath.exists(): # 如果配置文件不存在则创建配置文件
                self.cpath.write_text('{}')
            jstxt=self.cpath.read_text()
            try:
                self.conf:dict=json.loads(jstxt) # 如果文件格式损坏则重新生成空文件
                keys=self.conf.keys()
                for tip in self.titles.keys():
                    key=f'title.{tip}'
                    if key not in keys:continue
                    line=self.conf[key]
                    self.titles[tip].set_line(line)
            except:
                self.conf={}
                self.cpath.write_text('{}')
        
        self.search_title()


    def save_config(self): # 保存配置文件
        if not self.cache:return # 如果不生成缓存的话就直接返回吧
        with open(f"{self.cpath}",'w') as f:
            jstxt=json.dumps(self.conf,indent=4)
            f.write(jstxt)
    
    def search_title(self):
        """
        在文件中搜索需要的行号
        """
        from concurrent.futures import ThreadPoolExecutor
        from concurrent.futures._base import Future
        # 所有标题所在的行
        def sear_group(start:int): #每一个搜索的线程
            # print(f'search_group:{start}')
            rkeys=list(self.titles.keys()) # 总的key
            fkeys=[]
            for key in rkeys:
                title = self.titles[key]
                if title.line!=-1 and title.multi==False:continue # 如果已经有行数，而且不是多次匹配，则跳过
                fkeys.append(key)
            # keys=[key for key in keys if self.titles[key].line==-1] # 需要搜索的key,如果不是-1则代表不需要搜索
            if len(fkeys)==0:return
            for j in range(start,start+bsize):
                line=self.getline(j)
                if line=='':break
                for key in fkeys:
                    title:Title=self.titles[key]
                    if not title.judge(line):continue
                    self.titles[key].set_lnum(j)
                    # print(key,j)
                    if j>=self.titles[key].line:
                        # self.titles[key].line=j
                        self.conf[f'title_{key}']=j
                        self.save_config()
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
        lastLine='\n'.join(self.getlines(self.lineNum-10,self.lineNum))
        return 'Normal termination of Gaussian' in lastLine
    
    @lru_cache
    def get_atmXyzs(self)->np.ndarray:
        """原子坐标[n,3]"""
        result=self.read_coords()
        assert result is not None,"没有找到原子坐标"
        atmSyms,atmXyzs = result
        assert isinstance(atmXyzs,np.ndarray),'coord必须是np.ndarray类型'
        return atmXyzs

    @lru_cache
    def get_atmSyms(self)->list[str]:
        """原子符号[n]"""
        result=self.read_coords()
        assert result is not None,"没有找到原子符号"
        atmSyms,atmXyzs = result
        return atmSyms
    
    def get_nele(self)->tuple[int,int]: # 根据总核电荷数和电荷、自旋计算α、β电子数，不可行的
        lineNum=self.titles['nele'].line
        line=self.getline(lineNum)
        matchs=re.match(r'^ +(\d+) alpha electrons +(\d+) +beta electrons',line)
        assert matchs is not None,"没有找到电子数"
        nela,nelb=matchs.groups()
        return int(nela),int(nelb)
    
    def get_energy(self)->float:
        return self.read_energy()
    
    from pywfn.base.basis import BasisData,Basis
    def get_basis(self)->Basis:
        name=self.read_basisName()
        data=self.read_basisData()
        assert data is not None,"没有找到基组数据"
        basis=Basis()
        basis.name=name
        basis.data=data
        return basis
    
    
    def get_coefs(self)->Coefs:
        coefs=Coefs()
        atms,shls,syms,engs,occs,CM=self.read_CMs()
        coefs._atoAtms=atms
        coefs._atoShls=shls
        coefs._atoSyms=syms
        coefs._obtEngs=engs
        coefs._CM_raw=CM
        return coefs


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
                return symbols,np.array(coords)/0.529177 # 埃转为波尔
        raise AssertionError('没有读取到原子坐标')

    @lru_cache
    def read_multiy(self)->tuple[int,int]|None:
        """读取电荷和自选多重度"""
        res=re.findall(rf'Charge = +(-?\d) Multiplicity = (\d)',self.text)
        if res is None:
            printer.warn(f'{self.path} 没有读到电荷和自旋多重度!')
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
        读取基组数据，每个原子的都读出来
        原子类型
            电子层
        """
        from pywfn.base.basis import BasisData
        assert 'gfinput' in self.read_keyWrds(),'关键词应该包含：gfinput'
        titleNum=self.titles['basisData'].line
        symbols=self.get_atmSyms()
        if titleNum is None:return
        basisDatas:list[BasisData]=[]
        angDict={'S':0,'P':1,'D':2} #角动量对应的字典
        s1=rf'^ +(\d+) +\d+$'
        s2=r' ([SPD]+) +(\d+) \d.\d{2} +\d.\d{12}'
        s3=r'^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})'
        s4=' ****\n'
        elm=None
        shl=0 # 壳层
        angs=None
        exp=None
        coe=None
        # datas=[] # 存储基函数数据
        for i in range(titleNum+1,self.lineNum):
            line=self.getline(i)
            if re.search(s1,line) is not None:
                atm=re.search(s1,line).groups()[0] # type: ignore #第几个原子
                atm=int(atm)
                atmSym=symbols[atm-1] # 原子符号
                elem=elements[atmSym] # 原子类型
            elif re.search(s2,line) is not None:
                shellName,lineNum=re.search(s2,line).groups() # type: ignore
                angs=[angDict[s] for s in shellName] #角动量
                shl+=1
            elif re.search(s3,line) is not None:
                find=re.search(s3,line)
                assert find is not None,'正则匹配错误'
                numsStr=find.groups()[0]
                numStrs:list[str]=re.findall(r'-?\d.\d{10}D[+-]\d{2}',numsStr)
                nums:list[float]=[float(num.replace('D','E')) for num in numStrs]
                assert angs is not None,'angs is None'
                # assert atomic is not None,'atomic is None'
                if len(angs)==1:
                    alp,coe=nums
                    basisDatas.append(BasisData(atm,shl,angs[0],coe,alp))
                if len(angs)==2:
                    alp,coe1,coe2=nums
                    basisDatas.append(BasisData(atm,shl,angs[0],coe1,alp))
                    basisDatas.append(BasisData(atm,shl,angs[1],coe2,alp))
            elif line==s4 is not None:
                shl=0
                continue
            else:
                break
                # 排序一下
        basisDatas.sort(key=lambda b:(b.atm,b.shl,b.ang))
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
                    match=line[11:].strip()
                    shl,sym=match[0],match[1:] #这里壳层只能是整数哦
                    # print(line[11:],shl,sym)
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
        atms,shls,syms,engs,occs,CM
        """
        rdatas=[None]*6
        if self.cache:
            atms=self.load_fdata('atms.npy')
            shls=self.load_fdata('shls.npy')
            syms=self.load_fdata('angs.npy')
            engs=self.load_fdata('engs.npy')
            occs=self.load_fdata('occs.npy')
            CM  =self.load_fdata('CM.npy')
            rdatas=[atms,shls,syms,engs,occs,CM]

        if all(rdatas): # 所有读取的数据都不为空，有可能某些数据被用户删除
            lists:list[np.ndarray]=[atms,shls,syms,engs,occs] # type: ignore
            atms,shls,syms,engs,occs=[e.tolist() for e in lists]
            return atms,shls,syms,engs,occs,CM # type: ignore
        elif self.titles['coefs'].line!=-1:
            atms,shls,syms,engs,occs,CM=self.read_CM('coefs')
        elif self.titles['acoefs'].line!=-1:
            atmsA,shlsA,symsA,engsA,occsA,CMA=self.read_CM('acoefs')
            atmsB,shlsB,symsB,engsB,occsB,CMB=self.read_CM('bcoefs')
            atms=atmsA
            syms=symsA
            shls=shlsA
            engs=engsA+engsB
            occs=occsA+occsB
            CM=np.concatenate([CMA,CMB],axis=1)
        else:
            raise ValueError('没有找到轨道系数，请检查关键词是否完整')
        if self.cache:
            self.save_fdata('atms',atms)
            self.save_fdata('shls',shls)
            self.save_fdata('syms',syms)
            self.save_fdata('engs',engs)
            self.save_fdata('occs',occs)
            self.save_fdata('CM'  ,CM)
        # 将潜在的球谐的转为笛卡尔的形式
        # atmList,shlList,symList,CM=toCart(atms,shls,syms,CM)
        return atms,shls,syms,engs,occs,CM

    @lru_cache
    def read_SM(self):
        """读取重叠矩阵"""
        SM=self.read_Mat('overlap')
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
        res=re.search(r'Sum of Mulliken charges = +(-?\d.\d{5})',self.text).group(1) # type: ignore
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
            finds=re.findall(rf'-?\d+.\d+',line)
            freqs+=[float(e) for e in finds]
        return freqs
        
