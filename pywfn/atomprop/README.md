很多性质都是基于原子电荷的，目前可以计算的原子性质有以下：
- `原子电荷`有很多种，目前实现了`Mulliken`电荷和`lowdin`电荷，都在`atomCharge.py`中实现
- `自旋布局`是开壳层分子中原子的α电子数-β电子数，本质上是基于电荷计算的结果`atomSpin.py`
- `原子轨道能`是根据所有分子轨道中原子电荷的比例对分子轨道能量加权求和，结果依赖于电荷`atomEnergy.py`
- `fukui`函数和`parr`函数是不同价态分子`电荷`或`自旋`的差值，结果依赖于电荷`delProps.py`
- `轨道投影法`修改分子轨道系数矩阵，但还是要代入上述的基于电荷的分析方法`dirProps.py`
- `pi电子数`是利用轨道投影法指定方向为垂直于分子平面方向的方向电子数`piProps.py`
- `自由价`是根据与原子相邻的pi键级计算出的，与原子电荷没关系`freeValence.py`

## 原子轨道能
⚠以下公式都是我自己推导的
> 将每个分子轨道的电子分布到基函数上，每个基函数在不同轨道电子数比例×分子轨道能量=基函数轨道能，原子对应基函数的电子轨道能加和即为原子轨道能

密度矩阵`P`与重叠矩阵`S`的乘积的对角线元素为每个基函数的电子数量，以`b`代表基函数指标

$$N_b=(PS)_{b,b}$$

$$=P_{b,:}S_{:,b}$$

$$=\sum_{i=1}^{n_b} P_{b,i}S_{i,b}$$

$$
\sum_{i=1}^{n_{b}} \sum_{j=1}^{n_{o}}(C_{b,j}C_{j,b}) S_{i,b}
$$

其中`i`代表与该基函数`b`作用的基函数指标，`nb`代表基函数的数量，`j`代表分子轨道，`no`代表占据轨道数量，整个分子的电子数量为：

$$
N=\sum_{b=1}^{nb}\sum_{i=1}^{nb} \sum_{j=1}^{no}(C_{b,j}C_{j,b}) S_{i,b}
$$

第`b`个基函数在第`j`个分子轨道中的电子数量为：

$$
N_{b,j}=\sum_{i=1}^{nb}(C_{b,j}C_{j,b}) S_{i,b}
$$

其中`b`和`j`为常数；第`j`个轨道的电子数量为：

$$
N_{j}=\sum_{b=1}^{nb}\sum_{i=1}^{nb} (C_{b,j}C_{j,b}) S_{i,b}
$$

其中`j`为常数，该值等于1或2,记为`Ne`

第`j`个分子轨道中，第`b`个基函数电子数量的占比为：
$$
r_{b,j}=N_{b,j}/N_{e}
$$

则第`b`个基函数的电子在第`j`个分子轨道中的能量为：
$$
E_{b,j}=r_{b,j}\times e_j
$$

其中`ej`为第`j`个分子轨道的能量，则第`b`个基函数的电子的轨道能为：
$$
E_b=\sum_{j=1}^{no}E_{b,j}
$$

对原子`a`的所有基函数能量加和，得到原子的轨道能量

$$
E_a=\sum_{b∈a}E_b
$$

## hirshfeld电荷

在r点处的`promolecule density`为：
$$
\rho ^{pro}(r)=\sum_{i}\rho _i^{at}(r)
$$
其中$\rho _i^{at}$为基态原子密度。对于每一个原子，定义`sharing function`
$$
w_i(r)=\rho _i^{at}(r)/\rho ^{pro}(r)
$$
定义`bonded atom` i 的`charge density`：
$$
\rho _i^{b.a.}(r)=w_i(r)\rho ^{mol}(r)
$$
其中$\rho ^{mol}(r)$为`actual molecular density`

接下来有两种处理方式：
- 方式1：`atomic deformation density`：
$$
\delta \rho _i(r) =\rho _i^{b.a.}(r)-\rho _i^{at}(r)
$$
- 方式2：`molecule deformation density`
$$
\Delta  \rho (r) =\rho ^{mol}(r)-\rho ^{pro}(r)
$$
其中：
$$
\delta \rho_i (r) =w_i(r)\Delta \rho (r)
$$
可以计算原子电荷：
$$Q_i=-\int \rho _i^{b.i.}(r)dv$$

$$q_i=Q_i+Z_i$$

$$q_i=-\int \delta \rho _i(r)dv$$