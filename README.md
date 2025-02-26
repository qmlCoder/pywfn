# `pywfn` -- 基于python的波函数分析工具

文档： https://www.xiaofei911.top/mkdocs/pywfn/


## 依赖
```
numpy==2.1.1
rich==13.8.0
matplotlib==3.9.2
```
## 运行(CLI)
``` shell
python main.py
```

## 示例(API)
```python
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomProp import atomCharge

path="D:\BaiduSyncdisk\gfile\CnHn\C6H6.log" # 高斯输出文件的路径
reader=LogReader(path) # 实例化log文件读取器
mol=Mol(reader) # 实例化分子对象

caler=atomCharge.Calculator(mol) # 实例化原子电荷计算器，传入分子对象
result=caler.mulliken() # 计算mulliken电荷
print(result) # 打印结果
```

## 功能
![](./docs/pywfn_xmind.png)

gfortran -shared -ffree-form -ffree-line-length-none -fopenmp data.f90 march.f90 flib.f90 -o flib.dll

ifx -dll -free -qopenmp data.f90 march.f90 flib.f90 -o flib.dll