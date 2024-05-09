保存一些常用的数学公式

gto的计算量是很大的，所以可以尝试多种加速计算的方法
- 常规的numpy
- taichi
- julia
- f2py


加速计算做的尝试

使用numba的jit加速
使用fortran编译动态链接库
整个原子的gto作为高维矩阵一块算
原子轨道分层级计算
牺牲精度换速度
舍去不必要的点
系数非常小的也可不算

## python调用fortran
- 使用fortran的iso_c_binding
- 使用python的ctypes

iso_c_binding是fortran与c联用的工具

ctypes是python与fortran联用的工具

借用这两个工具，python就可以和fortran联用了

将fortran代码编译成动态链接库
使用python的ctypes调用动态链接库

> 关键问题是数据类型之间的对接，对c不太熟悉