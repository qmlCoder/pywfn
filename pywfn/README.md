这可能是一个长远的项目，以后计算别的东西也都可以放到这个项目里面  
所有计算部分都单独拎出来，保证程序在没有GUI的时候也可以当作脚本单独使用
# 项目结构 
## 基础类 bases
首先要有基础的类
- Mol 根据读取的文件生成分子对象
    - 其中atoms的索引是从1开始
    - 进行计算时主要使用的时Mol和Atom的方法和属性
    - 所有的计算都是基于分子类的，所以分子类应可以自定义
- Atom 根据读取的文件生成原子对象  
> 分子对象的信息可以自己指定的，也可以是从文件中读取的
## 读取器 readers
读取各种文件，如log/out或fchk
- LogReader 读取log/out文件
- FchReader 读取fch文件
## 计算器 calculates
实例化计算器，传入分子，计算的结果将保存到基础类的属性中  
获得属性 -> 计算指定任务 - 返回并设置属性
- 计算键级
    - 计算π键级(新方法) projectionPI
    - 计算π键级(老方法) selectionPI
    - 计算mayer键级 mayerOrder
    - 计算mayerπ键级 
- 原子属性
    - Mulliken电荷
## 工具类 tools
其它的一些对文件操作的函数
- Spliter 扫描/优化步骤的分割并单独文件导出
- Gif 指定信息，根据模板批量生成gjf文件
- SI 批量导出SI信息

## 通用函数 utils
有些计算需要使用大量的相同计算(通常为数学计算)，将这些计算保存在单独的文件`utils.py`中

# 获取actor对象
- 调用pv.add_mesh()时会返回actor对象，可以设置其不透明度

# 使用说明
```python
# 实例化一个Reader对象
reader=pywfn.Reader(path)
# 获取分子对象
mol=reader.mol
# 获取原子对象
atom=mol.atoms[1]

```
# 盘点一下程序中会出现的变量名以及缩写
> 变量名过长的话会导致程序冗长,而且函数会更长,所以尽量保持每一个含义长度在5以内  
总规则: 动词_形容词+名词
- 贡献 contribution cont
- 坐标 coordinates coord
- 系数 coefficients coeft
- 投影 projection projt
- 方向 direction 

一个原子的坐标为`coord`  
一个分子的所有原子的坐标为`coords`  
多个分子的原子坐标为`allCoords`
# 叨叨叨叨
每一个模块都有一个类，所以调用起来可能比较费劲：`pywfn.子包.子模块.类`  
一个模块只有一个类，感觉有点浪费，如果以后代码稳定了可以尝试把多个类写在同一个文件中  
由于想把每个功能都写在一个模块文件中，所以代码可能会显得比较冗杂,但是问题不大

## 支持的元素
> 有些元素根本用不着，所以不需要存储对不对

H Be C N O F P S Cl 

## 方向的确定

1埃=1.889726玻尔半径