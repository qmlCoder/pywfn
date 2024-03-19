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

from pywfn.data.basis import Basis
from pywfn.maths.gto import Gto
from pywfn.reader.lutils import Reader
from pywfn.utils import printer
from pywfn.data import elements

class Title:
    def __init__(self,line=None,pate='') -> None:
        self.line:int=line
        self.pate:str=pate
    
    def __repr__(self) -> str:
        return f'{self.line}: {self.pate}'

class LogReader(Reader):
    def __init__(self, path:str):
        Reader.__init__(self,path)
        if 'Normal termination of Gaussian' not in self.text:
            printer.wrong('文件未正常结束!')
        self.keyWords=re.search(r'^ # .+$',self.text,re.M).group()
        self.index=0 #从第一行向下搜索，搜索到需要的就停止
        self.titles={ #记录每一个title所在的行数
            'coords':Title(pate=r'(Input orientation|Standard orientation)'),
            'basis':Title(pate=r'Standard basis:'),
            'coefs':Title(pate=r'  Molecular Orbital Coefficients'),
            'acoefs':Title(pate=r'Alpha Molecular Orbital Coefficients'),
            'bcoefs':Title(pate=r'Beta Molecular Orbital Coefficients'),
            'overlap':Title(pate='\*\*\* Overlap \*\*\*'),
            'density':Title(pate=r'Density Matrix:'),
            'basisData':Title(pate=r'Overlap normalization'),
            'engs':Title(pate=r'Zero-point correction')
        }
        self.search_title()
        # self.read_basisData()
    
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
        atoms,layer,engs,type_,CM=self.read_CMs()
        return CM
    
    @lru_cache
    def get_obtAtoms(self) -> list[int]:
        atoms,layer,engs,type_,CM=self.read_CMs()
        return atoms
    
    @lru_cache
    def get_obtEngs(self) -> list[float]:
        atoms,layer,engs,type_,CM=self.read_CMs()
        return engs
    
    @lru_cache
    def get_obtLayer(self) -> list[str]:
        atoms,layer,engs,type_,CM=self.read_CMs()
        return layer
    
    def get_obtTypes(self) -> list[str]:
        atoms,layer,engs,type_,CM=self.read_CMs()
        return type_

    @lru_cache
    def get_SM(self)->np.ndarray:
        return self.read_SM()
    
    def get_charge(self) -> int:
        charge,spin=self.read_multiy()
        return charge
    
    def get_spin(self)->int:
        charge,spin=self.read_multiy()
        return spin
    
    def get_energy(self)->float:
        keys,nums=self.read_energy()
        return nums[-1]
    
    from pywfn.data import Basis
    def get_basis(self)->Basis:
        name=self.read_basisName()
        data=self.read_basisData()
        basis=Basis(name)
        basis.setData(data)
        return basis


    def search_title(self):
        """
        搜索到需要用到的标题行就停止
        """
        # 所有标题所在的行
        def sear_group(lines,start):
            for j,line in enumerate(lines):
                for key_ in self.titles.keys():
                    pate=self.titles[key_].pate
                    if re.search(pate,line) is None:continue #如果不匹配则跳过
                    self.titles[key_].line=start+j
                    printer.log(self.titles[key_])
        threads:list[threading.Thread]=[]
        for i in range(0,len(self.lines),100_000):
            lines=self.lines[i:i+100_000]
            t=threading.Thread(target=sear_group,args=(lines,i,))
            t.start()
            threads.append(t)
        for t in threads:
            t.join()

    @lru_cache
    def read_coords(self):
        '''读取原子坐标'''
        titleNum=self.titles['coords'].line
        if titleNum is None:
            printer.warn('没有读取到原子坐标')
            return
        s1=r' +\d+ +(\d+) +\d +(-?\d+.\d{6}) +(-?\d+.\d{6}) +(-?\d+.\d{6})'
        coords=[]
        symbols=[]
        for i in range(titleNum+5,len(self.lines)):
            line=self.lines[i]
            if re.search(s1, line) is not None:
                res=list(re.search(s1, line).groups())
                atomID=int(res[0])
                symbol=elements[atomID].symbol
                coord=[float(each) for each in res[1:]]
                symbols.append(symbol)
                coords.append(coord)
            else:
                return symbols,np.array(coords)

    @lru_cache
    def read_multiy(self):
        """读取电荷和自选多重度"""
        res=re.findall(f'Charge = +(-?\d) Multiplicity = (\d)',self.text)
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
            name=re.match(' Standard basis: (\S+) ',self.lines[titleNum]).group(1)
        return name
    
    @lru_cache
    def read_basisData(self):
        """
        读取基组数据
        原子类型
            电子层
        """
        assert 'gfinput' in self.keyWords,'关键词应该包含：gfinput'
        titleNum=self.titles['basisData'].line
        symbols=self.get_symbols()
        if titleNum is None:return
        basisData:dict[str:list[dict[str:list]]]={}
        ifRead=True
        angDict={'S':0,'P':1,'D':2} #角动量对应的字典
        s1='^ +(\d+) +\d+$'
        s2=r' ([SPD]+) +(\d+) \d.\d{2} +\d.\d{12}'
        s3=r'^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})'
        s4=' ****'
        for i in range(titleNum+1,len(self.lines)):
            line=self.lines[i]
            if re.search(s1,line) is not None:
                idx=re.search(s1,line).groups()[0]
                idx=int(idx)-1
                atomic=symbols[idx]
                if atomic not in basisData.keys():
                    basisData[atomic]=[] #shells
                    ifRead=True # 该元素是否已经读过
                else:
                    ifRead=False
            elif re.search(s2,line) is not None:
                if not ifRead:continue
                shellName,lineNum=re.search(s2,line).groups()
                ang=[angDict[s] for s in shellName] #角动量
                exp=[]
                # coe=[[]]*len(ang)
                coe=[[] for i in range(len(ang))]
                shell={'ang':ang,'exp':exp,'coe':coe}
                atomic:str=atomic
                basisData[atomic].append(shell)
            elif re.search(s3,line) is not None:
                if not ifRead:continue
                numsStr=re.search(s3,line).groups()[0]
                nums=re.findall(r'-?\d.\d{10}D[+-]\d{2}',numsStr)
                nums=[float(num.replace('D','E')) for num in nums]
                basisData[atomic][-1]['exp'].append(nums[0])
                basisData[atomic][-1]['coe'][0].append(nums[1])
                if len(ang)==2:
                    basisData[atomic][-1]['coe'][1].append(nums[2])
            elif line==s4 is not None:
                continue
            else:
                # self.mol.basis.setData(basisData)
                basisData=trans_basisData(basisData)
                return basisData
                break
        
    def read_summery(self):
        '''读取总结信息'''
        summerys=re.findall(r' \d[\\||]\d[\S\s]*?@',self.content)
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
    @lru_cache
    def read_CM_(self, title:str):  # 提取所有原子的轨道 自己写的代码自己看不懂真实一件可悲的事情,此函数逻辑复杂，要好好整明白
        s1='^(( +\d+){1,5}) *$'
        s2=r'^(( *(\(\w+\)--)?[OV]){1,5})$'
        s3=r'^ +Eigenvalues --(( +-?\d+.\d+){1,5})'
        s4='^ +\d+ +(\d+) +([A-Za-z]+) +(\d[A-Z]+)(( *-?\d+.\d+){1,5})$'
        # s5='^ +\d+ +(\d+[A-Za-z ]+)(( *-?\d+.\d+){1,5})$'
        s5='^ +\d+ +(\d+[A-Za-z]+ ?\+?-?\d?)(( *-?\d+.\d+){1,5})$'
        # titleNum=None
        # for i,line in enumerate(self.logLines):
        #     if title in line:
        #         titleNum=i
        titleNum=self.titles[title].line
        if titleNum is None:return None

        OrbitalAtom=[] # 系数矩阵每一行对应的原子
        OrbitalLayer=[] # 系数矩阵每一行对应的壳层符号
        OrbitalEngs=[] # 轨道本征值
        OrbitalType=[] # 轨道类型
        dataDict:dict[int:list[list[float]]]={}
        firstShow=True
        if titleNum is None:
            return
        for i in range(titleNum+1,len(self.lines)):
            line=self.lines[i]
            # print(line)
            if re.search(s1, line) is not None: #情况1
                pass
            elif re.search(s2, line) is not None: # 情况2，获得column
                OrbitalType += re.split(r' +', line.replace('\n',''))[1:] # 获取占据轨道还是非占据轨道
                
            elif re.search(s3, line) is not None: # 情况3，每个轨道本征值     
                line_data,_=re.search(s3,line).groups()
                line_data=re.findall('-?\d+.\d+', line_data)
                line_data=[float(each) for each in line_data]
                OrbitalEngs+=line_data
            elif re.search(s4,line) is not None: # 原子出现
                atomIDX,atomType,layer,line_data,_=re.search(s4,line).groups()
                atomIDX = int(atomIDX)
                nums=[float(each) for each in re.findall('-?\d+.\d+',line_data)]
                if atomIDX not in dataDict.keys():
                    dataDict[atomIDX]=[]
                else:
                    firstShow=False
                dataDict[atomIDX].append([])
                dataDict[atomIDX][-1].append(nums)
                if firstShow:OrbitalLayer.append(layer)
                if firstShow:OrbitalAtom.append(atomIDX)
                # self.OCdict[atomIDX].set(layer,nums)
            elif re.search(s5,line) is not None:
                layer,line_data,_=re.search(s5,line).groups()
                line_data=re.findall('-?\d+.\d+', line_data)
                nums=[float(each) for each in line_data]
                dataDict[atomIDX][-1].append(nums)
                if firstShow:OrbitalLayer.append(layer)
                if firstShow:OrbitalAtom.append(atomIDX)
                # self.OCdict[atomIDX].set(layer,nums)
            else: # 若不满足以上任意一种情况，说明已经查找完毕，则对收集到的数据进行处理
                # print(dataDict)
                printer.console.log(f'读取完成,i={i},line={line}')
                for atomic,matrics in dataDict.items():
                    for i,matrix in enumerate(matrics):
                        dataDict[atomic][i]=np.array(matrix)
                for atomic,matrics in dataDict.items():
                    dataDict[atomic]=np.concatenate(matrics,axis=1)
                CM=np.concatenate([m for m in dataDict.values()],axis=0)
                break
        printer.console.log(f'读取完成CM,shape={CM.shape}')
        assert CM.shape[0]==CM.shape[1],"CM需要为正方形矩阵"
        return OrbitalAtom,OrbitalLayer,OrbitalEngs,OrbitalType,CM
    
    @lru_cache
    def read_CM(self, title:str):  # 提取所有原子的轨道 自己写的代码自己看不懂真实一件可悲的事情,此函数逻辑复杂，要好好整明白
        assert 'pop=full' in self.keyWords,'关键词应包含：pop=full'
        find=re.search('NBasis *= *(\d+)',self.text).groups()[0]
        NBasis=int(find)
        NBlock=NBasis//5+(0 if NBasis%5==0 else 1)
        titleNum=self.titles[title].line
        
        if titleNum is None:return None
        blockLen=NBasis+3
        # print(titleNum,NBlock,blockLen)
        OrbitalAtom=[]
        OrbitalEngs=[]
        OrbitalType=[]
        OrbitalLaye=[]
        CM=np.zeros(shape=(NBasis,NBasis))
        for i,l in enumerate(range(titleNum+1,titleNum+NBlock*blockLen+1,blockLen)):
            # print(i,l,self.lines[l])
            types=self.lines[l+1][21:71]
            # print(types)
            types=[types[i:i+10].strip() for i in range(0,50,10)]
            types=[e for e in types if e!='']
            # print(types)
            # types=[float(e) for e in types]
            OrbitalType+=types
            engs =re.split(' +',self.lines[l+2][21:].strip())
            engs=[float(e) for e in engs]
            OrbitalEngs+=engs
            
            coefs=[line[21:] for line in self.lines[l+3:l+3+NBasis]]
            for j,line in enumerate(coefs):
                coef=re.findall('-?\d+.\d+',line)
                coef=[float(e) for e in coef]
                # print(i,j,coef)
                a=i*5
                b=i*5+len(coef)
                CM[j,i*5:i*5+len(coef)]=coef
            # print(coefs)

            if i==0:
                atomID=''
                for l2 in range(l+3,l+blockLen):
                    line=self.lines[l2][:15]
                    # print(l2,line)
                    obtLaye=line[12:].strip()
                    OrbitalLaye.append(obtLaye)
                    obtAtom=line[5:9].strip()
                    # print(obtAtom)
                    if obtAtom!='':
                        atomID=int(obtAtom)
                    OrbitalAtom.append(atomID)
                    # print(line,obtAtom,atomID)
        # print(NBasis)
        
        # print('OrbitalAtom',OrbitalAtom,len(OrbitalAtom))
        # print('OrbitalLaye',OrbitalLaye,len(OrbitalLaye))
        # print('OrbitalEngs',OrbitalEngs,len(OrbitalEngs))
        # print('OrbitalType',OrbitalType,len(OrbitalType))
        # print(CM.shape)
        return OrbitalAtom,OrbitalLaye,OrbitalEngs,OrbitalType,CM

    
    @lru_cache
    def read_CMs(self)->tuple[list,list,list,list,np.ndarray]:
        """系数矩阵"""
        if res:=self.read_CM('coefs'):
            atoms,layer,engs,type_,CM=res
        elif res:=self.read_CM('acoefs'):
            atomsA,layerA,engsA,type_A,CMA=res
            atomsB,layerB,engsB,type_B,CMB=self.read_CM('bcoefs')
            atoms=atomsA
            layer=layerA
            engs=engsA+engsB
            type_=type_A+type_B
            CM=np.concatenate([CMA,CMB],axis=1)
        return atoms,layer,engs,type_,CM

    @lru_cache
    def read_SM(self):
        """读取重叠矩阵"""
        s1='^( +\d+){1,5} *$'
        s2=' +\d+( +-?\d.\d{6}D[+-]\d{2}){1,5} *'
        lineDatas:dict[str:str]={}
        titleNum=self.titles['overlap'].line
        if titleNum is None:
            return None
        for i in range(titleNum+1,len(self.lines)):
            line=self.lines[i]
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
    
    @lru_cache
    def read_energy(self)->tuple[list[str],list[float]]:
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
        engDict={e:None for e in engList}
        lineNum=self.titles['engs'].line
        assert lineNum is not None,f"{self.path} 未读取到能量"
        for i in range(lineNum,len(self.lines)):
            line=self.lines[i]
            for each in engList:
                res=re.search(f'{each}=\s+(-?\d+\.\d+)',line)
                if res is not None:
                    re.findall('-?\d+\.\d+',line)
                    engDict[each]=res.groups()[0]
                    searhNum+=1
            if searhNum==len(engList):break
        engNums=[float(e) for e in engDict.values()]
        return engList,engNums

def trans_basisData(basisData:dict):
    """
    翻译基组数据
    原始的基组数据不好理解，将其转换为好理解的方式
    三重循环：原子、壳层、角动量
    每一个角动量对应二维数组[数所数量,2]
    2分别为收缩系数和高斯指数
    """
    newBasis=[]
    for atom,shells in basisData.items():
        for s,shell in enumerate(shells):
            angs:list[int]=shell['ang']
            exps:list[int]=shell['exp']
            coes:list[list[int]]=shell['coe']
            coes=np.array(coes)
            for a,ang in enumerate(angs):
                for e,exp in enumerate(exps):
                    coe=coes[a,e]
                    idx=elements[atom].idx
                    newBasis.append([idx,s,ang,exp,coe])
    return newBasis

def numlist(l):
    """将列表字符串转为数字"""
    return [float(e) for e in l]