# `pywfn` -- 基于python的波函数分析工具

文档： https://www.xiaofei911.top/mkdocs/pywfn/


## 依赖
```
numpy>=2.1.1
```
> python >= 3.10

## 运行(CLI)
``` shell
python main.py
```

## 示例(API)
```python
from pywfn.base import Mole
from pywfn.props.atomprop import charge

path="C6H6.log" # 高斯输出文件的路径
mole=Mole.from_file(path) # 实例化分子对象

caler=charge.Calculator(mole) # 实例化原子电荷计算器，传入分子对象
result=caler.mulliken() # 计算mulliken电荷
print(result) # 打印结果
```

## 功能
![](./docs/pywfn_xmind.png)
