use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
#[derive(Clone)]
pub struct Ato {
    pub core: rswfn::base::Ato,
}

#[pymethods]
impl Ato {
    /// 创建一个原子轨道
    ///
    /// **参数**
    /// - `atm`: 所属原子索引
    /// - `shl`: 所属壳层索引
    /// - `ang`: 角动量量子数
    /// - `sym`: 轨道符号
    /// - `loc`: 起始位置
    /// - `form`: 基函数形式
    /// - `alps`: 高斯指数列表
    /// - `coes`: 收缩系数列表
    #[new]
    pub fn new(
        atm: usize,
        shl: usize,
        ang: usize,
        sym: String,
        loc: usize,
        form: String,
        alps: Vec<f64>,
        coes: Vec<f64>,
    ) -> Self {
        Self {
            core: rswfn::base::Ato::new(atm, shl, ang, sym, loc, form, &alps, &coes),
        }
    }

    /// 所属原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }

    /// 所属壳层索引
    #[getter]
    pub fn shl(&self) -> usize {
        self.core.shl
    }

    /// 角动量量子数
    #[getter]
    pub fn ang(&self) -> usize {
        self.core.ang
    }

    /// 轨道符号
    #[getter]
    pub fn sym(&self) -> String {
        self.core.sym.clone()
    }

    /// 起始位置
    #[getter]
    pub fn loc(&self) -> usize {
        self.core.loc
    }

    /// 基函数形式
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }

    /// 高斯指数列表
    #[getter]
    pub fn alps(&self) -> Vec<f64> {
        self.core.alps.clone()
    }

    /// 收缩系数列表
    #[getter]
    pub fn coes(&self) -> Vec<f64> {
        self.core.coes.clone()
    }

    /// 获取归一化后的收缩系数
    pub fn coes_norm(&self) -> Vec<f64> {
        self.core.coes_norm().clone()
    }

    /// 获取量子数 (l, m, n)
    pub fn lmn(&self) -> [usize; 3] {
        *self.core.lmn()
    }
}

#[pyclass]
pub struct Atos {
    pub core: rswfn::base::Atos,
}

#[pymethods]
impl Atos {
    /// 返回指定基函数形式的原子轨道数量
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    pub fn len(&self, form: &str) -> usize {
        self.core.len(form)
    }

    /// 获取指定基函数形式下各壳层对应的原子索引列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn atms(&self, form: &str) -> Vec<usize> {
        self.core.atms(form)
    }

    /// 获取指定基函数形式下各壳层对应的壳层索引列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn shls(&self, form: &str) -> Vec<usize> {
        self.core.shls(form)
    }

    /// 获取指定基函数形式下各壳层的符号表示列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn syms(&self, form: &str) -> Vec<String> {
        self.core.syms(form)
    }

    /// 获取指定基函数形式下各轨道的起始位置列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn locs(&self, form: &str) -> Vec<usize> {
        self.core.locs(form)
    }
}

#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct Shl {
    pub core: rswfn::base::Shl,
}

#[pymethods]
impl Shl {
    /// 创建一个壳层
    ///
    /// **参数**
    /// - `atm`: 所属原子索引
    /// - `shl`: 壳层索引
    /// - `ang`: 角动量量子数
    /// - `form`: 基函数形式
    /// - `alps`: 高斯指数列表
    /// - `coes`: 收缩系数列表
    #[new]
    pub fn new(
        atm: usize,
        shl: usize,
        ang: usize,
        form: String,
        alps: Vec<f64>,
        coes: Vec<f64>,
    ) -> Self {
        Self {
            core: rswfn::base::Shl::new(atm, shl, ang, form, alps, coes),
        }
    }

    /// 获取该壳层所属的原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }
    /// 获取该壳层的索引编号
    #[getter]
    pub fn shl(&self) -> usize {
        self.core.shl
    }
    /// 获取该壳层的角动量量子数
    #[getter]
    pub fn ang(&self) -> usize {
        self.core.ang
    }
    /// 获取该壳层的起始轨道索引位置
    #[getter]
    pub fn loc(&self) -> usize {
        self.core.loc()
    }
    /// 获取该壳层的基函数形式
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }
    /// 获取该壳层的高斯指数（alpha 值）列表
    #[getter]
    pub fn alps(&self) -> Vec<f64> {
        self.core.alps.clone()
    }
    /// 获取该壳层的收缩系数列表
    #[getter]
    pub fn coes(&self) -> Vec<f64> {
        self.core.coes.clone()
    }

    /// 获取归一化后的收缩系数列表
    pub fn coes_norm(&self) -> Vec<f64> {
        self.core.coes_norm()
    }
}

#[pyclass]
pub struct Shls {
    pub core: rswfn::base::Shls,
}

#[pymethods]
impl Shls {
    /// 创建一个壳层集合
    ///
    /// **参数**
    /// - `shls`: 壳层列表
    #[new]
    pub fn new(shls: Vec<Shl>) -> Self {
        let shls = shls.into_iter().map(|shl| shl.core).collect();
        Self {
            core: rswfn::base::Shls::new(shls),
        }
    }

    /// 获取指定索引的壳层
    ///
    /// **参数**
    /// - `i`: 壳层索引
    ///
    /// **返回**
    /// - 壳层对象
    pub fn get(&self, i: usize) -> Shl {
        Shl {
            core: self.core.get(i).clone(),
        }
    }

    /// 返回壳层总数
    pub fn len(&self) -> usize {
        self.core.len()
    }

    /// 获取指定基函数形式的轨道上三角索引对列表
    ///
    /// 返回每对轨道的 (i, j) 索引，用于构建重叠矩阵等
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.core.uplo(form).clone()
    }

    /// 获取指定基函数形式的基函数数量
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn nato(&self, form: &str) -> usize {
        self.core.nato(form)
    }

    /// 获取指定基函数形式的原子对索引列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.core.atm_uplo(form)
    }
}

/// 单个高斯基函数（primitive）的原始数据
#[pyclass]
#[derive(Clone)]
pub struct Bas {
    pub core: rswfn::base::Bas,
}

#[pymethods]
impl Bas {
    /// 创建一个高斯基函数
    ///
    /// **参数**
    /// - `atm`: 所属原子索引
    /// - `shl`: 所属壳层索引
    /// - `ang`: 角动量量子数
    /// - `alp`: 高斯指数
    /// - `coe`: 收缩系数
    #[new]
    pub fn new(atm: usize, shl: usize, ang: usize, alp: f64, coe: f64) -> Self {
        Self {
            core: rswfn::base::Bas::new(atm, shl, ang, alp, coe),
        }
    }

    /// 所属原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }

    /// 所属壳层索引
    #[getter]
    pub fn shl(&self) -> usize {
        self.core.shl
    }

    /// 角动量量子数
    #[getter]
    pub fn ang(&self) -> usize {
        self.core.ang
    }

    /// 高斯指数
    #[getter]
    pub fn alp(&self) -> f64 {
        self.core.alp
    }

    /// 收缩系数
    #[getter]
    pub fn coe(&self) -> f64 {
        self.core.coe
    }
}

#[derive(Clone)]
#[pyclass(from_py_object)]
pub struct Basis {
    pub core: rswfn::base::Basis,
}

#[pymethods]
impl Basis {
    /// 创建一个基组
    ///
    /// **参数**
    /// - `shls`: 壳层列表
    /// - `name`: 基组名称，可选
    #[new]
    #[pyo3(signature = (shls=None, name=None))]
    pub fn new(shls: Option<Vec<Shl>>, name: Option<String>) -> Self {
        let shls = shls
            .unwrap_or_default()
            .into_iter()
            .map(|shl| shl.core)
            .collect();
        Self {
            core: rswfn::base::Basis::new(shls, name),
        }
    }

    /// 重建基组
    ///
    /// **参数**
    /// - `shls`: 壳层列表
    pub fn build(&mut self, shls: Vec<Shl>) {
        let shls = shls.into_iter().map(|shl| shl.core).collect();
        self.core.build(shls)
    }

    /// 判断基组是否为空
    pub fn is_none(&self) -> bool {
        self.core.is_none()
    }

    /// 获取积分使用的基函数形式
    pub fn int_form(&self) -> String {
        self.core.int_form().to_string()
    }

    /// 把基组统一到指定基函数形式
    ///
    /// **参数**
    /// - `form`: 目标基函数形式
    pub fn unify_form(&self, form: &str) -> Basis {
        Basis {
            core: self.core.unify_form(form),
        }
    }

    /// 获取指定原子的所有壳层
    ///
    /// **参数**
    /// - `atm`: 原子索引
    pub fn get_by_atm(&self, atm: usize) -> Vec<Shl> {
        self.core
            .get_by_atm(atm)
            .into_iter()
            .map(|core| Shl { core })
            .collect()
    }

    /// 获取所有原子轨道的键
    pub fn ato_keys(&self) -> Vec<String> {
        self.core.ato_keys()
    }

    /// 获取基函数数量
    pub fn nbas(&self) -> usize {
        self.core.nbas()
    }

    /// 获取壳层数量
    pub fn nshl(&self) -> usize {
        self.core.nshl()
    }

    /// 判断两个基组是否相同
    ///
    /// **参数**
    /// - `other`: 另一个基组
    pub fn equal(&self, other: Basis) -> bool {
        self.core.equal(&other.core)
    }

    /// 返回基组的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }

    /// 获取原子轨道壳层索引信息
    pub fn atos(&self) -> Atos {
        Atos {
            core: self.core.atos().clone(),
        }
    }

    /// 获取所有壳层对象
    pub fn shls(&self) -> Shls {
        Shls {
            core: self.core.shls.clone(),
        }
    }

    /// 获取指定基函数形式的重叠矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式
    /// - `mole`: 分子对象
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_smat(&self, py: Python, form: &str, mole: Mole) -> Py<PyArray2<f64>> {
        let smat = self.core.get_smat(form, &mole.core).to_owned();
        smat.into_pyarray(py).unbind()
    }

    /// 获取指定基函数形式的原子对索引列表
    ///
    /// 返回每对相互作用的原子索引 (i, j)，用于构建原子基的重叠矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        let uplo = self.core.atm_uplo(form).to_owned();
        uplo
    }
}
