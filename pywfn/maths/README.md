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

```shell
gfortran -shared flib.f90 -o flib.dll
```
- 使用fortran的iso_c_binding
- 使用python的ctypes

iso_c_binding是fortran与c联用的工具

ctypes是python与fortran联用的工具

借用这两个工具，python就可以和fortran联用了

将fortran代码编译成动态链接库
使用python的ctypes调用动态链接库

> 关键问题是数据类型之间的对接，对c不太熟悉

一般的阐述只有三种
- int,c_int
- float,c_double
- np.ndarray,POINTER(c_double)

- 传入的数据全部要使用np.float64类型
- 当python传入的参数数量与Fortran中的对应不上，python会直接停止
- fortran代码执行一部分就停止了
- 对函数绑定C命名时需要全为小写字母
- 不要单独的写函数，而是要写到模块里
- c_long:int32,c_float:float32
- 中间变量在函数中要先初始化为0
- 使用高级索引后的数组作为参数时需要copy一下

## 波函数的计算
高斯基函数gtf用$γ$表示，表达式为：
$$
γ=N.x^l.y^m.z^n.e^{-αr^2}
$$
其中：$r^2=x^2+y^2+z^2$，$α$为gtf的指数，N为归一化系数：
$$
N=\left ( \frac{2\alpha }{\pi }  \right ) ^{3/4} \sqrt{\frac{(4\alpha )^{l+m+n}}{(2l-1)!!(2m-1)!!(2n-1)!!} } 
$$

收缩基函数（原子轨道）cgf用$φ$表示，表达式为：
$$
\psi =\sum_{i=1}^{nc} c_i\gamma _i
$$
其中$c$为收缩系数，$nc$为收缩次数

分子轨道为原子轨道的线性组合，表达式为：
$$
\Psi=\sum_{i=1}^{nb} C_i\psi _i
$$
其中$nb$为基函数数量，$C$为轨道系数

## 数值重叠矩阵
重叠矩阵矩阵元表达式为：
$$
S_{u,v}=\int \psi _u\psi _v dr
$$
对空间采点，使用数值的方式为：
$$
S_{u,v}\sum_{i=1}^{np} \psi _u(r_i)\psi _v(r_i)W_i
$$
其中$np$为采点的数量，$W_i$为该点的权重