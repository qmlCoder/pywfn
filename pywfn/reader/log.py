"""
此脚本用来提取高斯输出文件的信息
高斯的输出文件包含迭代信息(结构优化、扫描等)
但是封装成分子对象之后就只有一个信息了
所以迭代信息只能在reader对象中存在,且不是默认属性
输出文件中存储的系数可能是球谐的，但是要转为笛卡尔的才方便使用
"""


from pywfn.base.basis import Basis
from pywfn.base.coefs import Coefs
from pywfn.base.geome import Geome
from pywfn.utils import printer
from pywfn.data.elements import elements
from pywfn import reader
from pywfn import core
import textwrap


class LogReader(reader.Reader):
    def __init__(self, path:str,cache:bool=False) -> None:
        super().__init__(path,cache)
        self.reader=core.reader.LogReader(path) # type: ignore
    
    # @property
    # def normalEnd(self)->bool:
    #     lastLine='\n'.join(self.getlines(self.lineNum-10,self.lineNum))
    #     return 'Normal termination of Gaussian' in lastLine
    
    
    def get_geome(self)->"Geome":
        """获取分子几何信息"""
        
        geome_core=self.reader.get_geome()
        geome=Geome()
        geome.core=geome_core
        return geome
    
    def get_basis(self)->"Basis":
        basis_core=self.reader.get_basis()
        basis=Basis()
        basis.core=basis_core
        return basis
    
    
    def get_coefs(self)->"Coefs":
        coefs_core=self.reader.get_coefs()
        coefs=Coefs()
        coefs.core=coefs_core
        return coefs
    
    # def get_nele(self)->tuple[int,int]: # 根据总核电荷数和电荷、自旋计算α、β电子数，不可行的
    #     lineNum=self.titles['nele'].line
    #     line=self.getline(lineNum)
    #     matchs=re.match(r'^ +(\d+) alpha electrons +(\d+) +beta electrons',line)
    #     assert matchs is not None,"没有找到电子数"
    #     nela,nelb=matchs.groups()
    #     return int(nela),int(nelb)
    
    # def get_energy(self)->float:
    #     return self.read_energy()
    
    # from pywfn.base.basis import BasisData,Basis
    


    # def read_keyWrds(self):
    #     """读取关键字"""
    #     titleNum=self.titles['keyWards'].line
    #     keyWards=''
    #     for i in range(titleNum,titleNum+3):
    #         line=self.getline(i,False).strip()
    #         if line=='-'*70:break
    #         keyWards+=line
    #     return keyWards

    # @lru_cache
    # def read_coords(self)->tuple[list[str],np.ndarray]|None:
    #     '''读取原子坐标'''
        
    #     titleNum=self.titles['coords'].line
    #     assert titleNum!=-1,"没有坐标对应的行"
    #     if titleNum is None:
    #         printer.warn('没有读取到原子坐标')
    #         return None
    #     s1=r' +\d+ +(\d+) +\d +(-?\d+.\d{6}) +(-?\d+.\d{6}) +(-?\d+.\d{6})'
    #     coords:list[list[float]]=[]
    #     symbols=[]
    #     for i in range(titleNum+5,self.lineNum):
    #         line=self.getline(i)
    #         find=re.search(s1, line)
    #         if find is not None:
    #             res=list(find.groups())
    #             atomID=int(res[0])
    #             symbol=elements[atomID].symbol
    #             coord=[float(each) for each in res[1:]]
    #             symbols.append(symbol)
    #             coords.append(coord)
    #         else:
    #             return symbols,np.array(coords)/0.529177 # 埃转为波尔
    #     raise AssertionError('没有读取到原子坐标')

    # @lru_cache
    # def read_multiy(self)->tuple[int,int]|None:
    #     """读取电荷和自选多重度"""
    #     res=re.findall(rf'Charge = +(-?\d) Multiplicity = (\d)',self.text)
    #     if res is None:
    #         printer.warn(f'{self.path} 没有读到电荷和自旋多重度!')
    #     else:
    #         charge,multiy=res[0]
    #         return int(charge),int(multiy)

    # @lru_cache 
    # def read_basisName(self)->str:
    #     """读取基组名"""
    #     titleNum=self.titles['basis'].line
    #     if titleNum is None:
    #         printer.warn('未读取到基组名')
    #         name='unll'
    #     else:
    #         find=re.match(rf' Standard basis: (\S+) ',self.getline(titleNum))
    #         if find is None:
    #             name='null'
    #         else:
    #             name=find.group(1)
    #     return name
    
    # @lru_cache
    # def read_basisData(self):
    #     """
    #     读取基组数据，每个原子的都读出来
    #     原子类型
    #         电子层
    #     """
    #     from pywfn.base.basis import BasisData
    #     assert 'gfinput' in self.read_keyWrds(),'关键词应该包含：gfinput'
    #     titleNum=self.titles['basisData'].line
    #     corods=self.read_coords()
    #     assert corods is not None,'没有原子坐标'
    #     symbols,_=corods
    #     if titleNum is None:return
    #     basisDatas:list[BasisData]=[]
    #     angDict={'S':0,'P':1,'D':2} #角动量对应的字典
    #     s1=rf'^ +(\d+) +\d+$'
    #     s2=r' ([SPD]+) +(\d+) \d.\d{2} +\d.\d{12}'
    #     s3=r'^ +(( +-?\d.\d{10}D[+-]\d{2}){2,3})'
    #     s4=' ****\n'
    #     elm=None
    #     shl=0 # 壳层
    #     angs=None
    #     exp=None
    #     coe=None
    #     # datas=[] # 存储基函数数据
    #     for i in range(titleNum+1,self.lineNum):
    #         line=self.getline(i)
    #         if re.search(s1,line) is not None:
    #             atm=re.search(s1,line).groups()[0] # type: ignore #第几个原子
    #             atm=int(atm)
    #             atmSym=symbols[atm-1] # 原子符号
    #             elem=elements[atmSym] # 原子类型
    #         elif re.search(s2,line) is not None:
    #             shellName,lineNum=re.search(s2,line).groups() # type: ignore
    #             angs=[angDict[s] for s in shellName] #角动量
    #             shl+=1
    #         elif re.search(s3,line) is not None:
    #             find=re.search(s3,line)
    #             assert find is not None,'正则匹配错误'
    #             numsStr=find.groups()[0]
    #             numStrs:list[str]=re.findall(r'-?\d.\d{10}D[+-]\d{2}',numsStr)
    #             nums:list[float]=[float(num.replace('D','E')) for num in numStrs]
    #             assert angs is not None,'angs is None'
    #             # assert atomic is not None,'atomic is None'
    #             if len(angs)==1:
    #                 alp,coe=nums
    #                 basisDatas.append(BasisData(atm,shl,angs[0],alp,coe))
    #             if len(angs)==2:
    #                 alp,coe1,coe2=nums
    #                 basisDatas.append(BasisData(atm,shl,angs[0],alp,coe1))
    #                 basisDatas.append(BasisData(atm,shl,angs[1],alp,coe2))
    #         elif line==s4 is not None:
    #             shl=0
    #             continue
    #         else:
    #             break
    #             # 排序一下
    #     basisDatas.sort(key=lambda b:(b.atm,b.shl,b.ang))
    #     return basisDatas
        
    # def read_summery(self):
    #     '''读取总结信息'''
    #     summerys=re.findall(r' \d[\\||]\d[\S\s]*?@',self.text)
    #     summery=summerys[-1]
    #     summery=''.join([each[1:] for each in summery.split('\n')])
    #     summeryList=summery.replace('\n','').replace('\\\\','||').split('||')
    #     Infos,KeyWords,JobTitle,Coords=summeryList[:4]
    #     _,_,_,jobType,method,basisSet,molFormula,user,date,_=Infos.replace('\\','|').split('|')
    #     summery={
    #         'basisSet':basisSet,
    #         'Coords':np.array([[float(num.replace(' ','')) for num in each.split(',')[1:]] for each in Coords.replace('\\','|').split('|')[1:]])
    #     }
    #     return summery
    
    
    # @lru_cache
    # def read_CM(self, title:str)->tuple[list[int],list[int],list[str],list[float],list[bool],np.ndarray]:  # 提取所有原子的轨道 自己写的代码自己看不懂真实一件可悲的事情,此函数逻辑复杂，要好好整明白
        
    #     keyWards=self.read_keyWrds()
    #     assert 'pop=full' in keyWards,'关键词应包含：pop=full'
    #     find=re.search(r'NBasis *= *(\d+)',self.text)
    #     assert find is not None,"没有找到基函数数量!!"
    #     NBasis=int(find.groups()[0])
    #     NBlock=NBasis//5+(0 if NBasis%5==0 else 1)
    #     titleNum=self.titles[title].line
        
    #     assert titleNum!=-1,"没有轨道系数!!"
    #     blockLen=NBasis+3 #一块数据行数，在log文件的输出中，轨道系数是按照每一块五列来分块的
    #     # print(titleNum,NBlock,blockLen)
    #     ObtAtms=[] # 原子轨道对应的原子
    #     ObtShls=[] # 占据原子的第多少层轨道
    #     ObtSyms=[] # 原子轨道对应的角动量
    #     ObtEngs=[] # 分子轨道能量
    #     ObtOccs=[] # 占据/非占据
        
    #     CM=np.zeros(shape=(NBasis,NBasis))
    #     for i,l in enumerate(range(titleNum+1,titleNum+NBlock*blockLen+1,blockLen)):
    #         ltex=self.getline(l+1)
    #         occs=self.getline(l+1)[21:71] # 占据类型的起止位置
    #         occs=[occs[i:i+10].strip() for i in range(0,50,10)]
    #         occs=[e for e in occs if e!='']
    #         occs=[True if e[-1]=='O' else False for e in occs]
    #         ObtOccs+=occs
    #         # engs =re.split(' +',self.getline(l+2)[21:].strip())
    #         engs =textwrap.wrap(self.getline(l+2)[21:].strip(),10)
    #         line=self.getline(l+2)
    #         engs=[float(e) for e in engs]
    #         ObtEngs+=engs
            
    #         coefs=[line[21:] for line in self.getlines(l+3,l+3+NBasis)]
    #         for j,line in enumerate(coefs):
    #             coef=re.findall(r'-?\d+.\d+',line)
    #             coef=[float(e) for e in coef]
    #             CM[j,i*5:i*5+len(coef)]=coef

    #         if i==0: #只在第一块记录行信息
    #             atomID=''
    #             for l2 in range(l+3,l+blockLen):
    #                 line=self.getline(l2)[:16] # 行信息(角动量，对应原子)的起止位置
    #                 match=line[11:].strip()
    #                 shl,sym=match[0],match[1:] #这里壳层只能是整数哦
    #                 # print(line[11:],shl,sym)
    #                 ObtShls.append(int(shl))
    #                 ObtSyms.append(sym)
    #                 obtAtom=line[5:9].strip()
    #                 if obtAtom!='':
    #                     atomID=int(obtAtom) # 更改当前行的原子
    #                 ObtAtms.append(atomID)
        
    #     return ObtAtms,ObtShls,ObtSyms,ObtEngs,ObtOccs,CM

    # @lru_cache
    # def read_CMs(self)->tuple[list,list,list,list,list,np.ndarray]:
    #     """
    #     获取轨道系数及相关信息
    #     atms,shls,syms,engs,occs,CM
    #     """
    #     rdatas=[None]*6
    #     if self.cache:
    #         atms=self.load_fdata('atms.npy')
    #         shls=self.load_fdata('shls.npy')
    #         syms=self.load_fdata('angs.npy')
    #         engs=self.load_fdata('engs.npy')
    #         occs=self.load_fdata('occs.npy')
    #         CM  =self.load_fdata('CM.npy')
    #         rdatas=[atms,shls,syms,engs,occs,CM]

    #     if all(rdatas): # 所有读取的数据都不为空，有可能某些数据被用户删除
    #         lists:list[np.ndarray]=[atms,shls,syms,engs,occs] # type: ignore
    #         atms,shls,syms,engs,occs=[e.tolist() for e in lists]
    #         return atms,shls,syms,engs,occs,CM # type: ignore
    #     elif self.titles['coefs'].line!=-1:
    #         atms,shls,syms,engs,occs,CM=self.read_CM('coefs')
    #     elif self.titles['acoefs'].line!=-1:
    #         atmsA,shlsA,symsA,engsA,occsA,CMA=self.read_CM('acoefs')
    #         atmsB,shlsB,symsB,engsB,occsB,CMB=self.read_CM('bcoefs')
    #         atms=atmsA
    #         syms=symsA
    #         shls=shlsA
    #         engs=engsA+engsB
    #         occs=occsA+occsB
    #         CM=np.concatenate([CMA,CMB],axis=1)
    #     else:
    #         raise ValueError('没有找到轨道系数，请检查关键词是否完整')
    #     if self.cache:
    #         self.save_fdata('atms',atms)
    #         self.save_fdata('shls',shls)
    #         self.save_fdata('syms',syms)
    #         self.save_fdata('engs',engs)
    #         self.save_fdata('occs',occs)
    #         self.save_fdata('CM'  ,CM)
    #     # 将潜在的球谐的转为笛卡尔的形式
    #     # atmList,shlList,symList,CM=toCart(atms,shls,syms,CM)
    #     return atms,shls,syms,engs,occs,CM

    # @lru_cache
    # def read_SM(self):
    #     """读取重叠矩阵"""
    #     SM=self.read_Mat('overlap')
    #     return SM

    # def read_Mat(self,title:str)->np.ndarray:
    #     """读取矩阵"""
    #     s1=r'^( +\d+){1,5} *$'
    #     s2=r' +\d+( +-?\d.\d{6}D[+-]\d{2}){1,5} *'
    #     lineDatas:dict[str,str]={}
    #     titleNum=self.titles[title].line
    #     for i in range(titleNum+1,self.lineNum):
    #         line=self.getline(i)
    #         if re.match(s1, line) is not None:
    #             pass
    #         elif re.match(s2, line) is not None:
    #             idx=line[:8]
    #             if idx not in lineDatas.keys():
    #                 lineDatas[idx]=''
    #             nums=line[8:].replace('D','e')
    #             lineDatas[idx]+=f'{nums} '
    #         else:
    #             break
    #     size=len(lineDatas.keys())
    #     SM_=np.zeros((size,size))
    #     for idx,each in lineDatas.items():
    #         each:str=each
    #         nums=re.split(' +',each.strip())
    #         nums=[float(num) for num in nums]
    #         SM_[int(idx)-1,:len(nums)]=nums
    #     return np.tril(SM_)+np.tril(SM_,-1).T

    # @lru_cache
    # def read_charge(self)->float:
    #     res=re.search(r'Sum of Mulliken charges = +(-?\d.\d{5})',self.text).group(1) # type: ignore
    #     return float(res)
    
    # def read_energy(self)->float:
    #     engs=re.findall(r'SCF Done: +E\(.*\) += +(-?\d+.\d+)',self.text)
    #     if engs:return float(engs[-1])
    #     printer.warn('未读取到能量，使用0代替')
    #     return 0.0

    # @lru_cache
    # def read_energys(self)->list[float]:
    #     """返回表格中对应的这8个能量，只返回数值"""
    #     engList = [
    #         'Zero-point correction', 
    #         'Thermal correction to Energy', 
    #         'Thermal correction to Enthalpy',
    #         'Thermal correction to Gibbs Free Energy', 
    #         'Sum of electronic and zero-point Energies',
    #         'Sum of electronic and thermal Energies', 
    #         'Sum of electronic and thermal Enthalpies',
    #         'Sum of electronic and thermal Free Energies'
    #         ]
    #     searhNum=0
    #     engDict={e:0.0 for e in engList}
    #     lineNum=self.titles['engs'].line
    #     if lineNum is None:printer.log("{self.path} 未读取到能量")
    #     for i in range(lineNum,self.lineNum):
    #         line=self.getline(i)
    #         for each in engList:
    #             res=re.search(rf'{each}=\s+(-?\d+\.\d+)',line)
    #             if res is not None:
    #                 re.findall(r'-?\d+\.\d+',line)
    #                 engDict[each]=float(res.groups()[0])
    #                 searhNum+=1
    #         if searhNum==len(engList):break
    #     engNums=[float(e) for e in engDict.values()]
    #     return engNums
    
    # def read_freqs(self):
    #     """读取频率"""
    #     freqs=[]
    #     find = re.findall(r'^\s+Frequencies\s--\s+(-?\d+\.\d+.+)$', self.text, flags=re.M)
    #     if find is None:return freqs
    #     for line in find:
    #         finds=re.findall(rf'-?\d+.\d+',line)
    #         freqs+=[float(e) for e in finds]
    #     return freqs
        
