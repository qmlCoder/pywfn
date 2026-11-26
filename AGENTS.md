这是rswfn的python包装，本质上是用python给我写的rust库套个壳
技术上使用 pyo3+maturin 架构

src_rs是先用pyo3写的python库，会生成一个_core.pyd文件，但是由于没有类型及代码提示因此需要再加一层python封装
src_py/python是最终的python代码，内部的所有逻辑都通过pywfn._core实现，函数的参数以及返回值都有应该有类型标注

src_rs不是直接使用的，因此代码不需要有注释
src_py/python是最终需要用到的，因此要有注释，注释来源参考rswfn

rswfn的位置可以在cargo.toml中找到

src_rs以及src_py/python内的结构应该与rswfn的结构保持一致

tests是测试代码，每次修改后要保证所有测试通过

python库通过maturin构建，`maturin develop`是开发版，测试时使用这个，测试通过之后使用`maturin build --release`生成发布版

发布时，更新一下版本号，发布到`pypi`以及`github release`中
