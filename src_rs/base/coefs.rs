use pyo3::prelude::*;

use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use rswfn;

use crate::{base::Mole, utils};

#[derive(Clone)]
#[pyclass(from_py_object)]
pub struct Coefs {
    pub core: rswfn::base::Coefs,
}

#[pymethods]
impl Coefs {
    /// 创建一个轨道系数对象
    ///
    /// **参数**
    /// - `ato_atms`: 每个基函数所属的原子索引
    /// - `ato_shls`: 每个基函数所属的壳层索引
    /// - `ato_syms`: 每个基函数的符号
    /// - `obt_engs`: 轨道能量列表
    /// - `obt_occs`: 轨道占据数列表
    /// - `cmat`: 轨道系数矩阵（NumPy 二维数组）
    #[new]
    pub fn new(
        ato_atms: Vec<usize>,
        ato_shls: Vec<usize>,
        ato_syms: Vec<String>,
        obt_engs: Vec<f64>,
        obt_occs: Vec<f64>,
        cmat: Bound<'_, PyArray2<f64>>,
    ) -> Self {
        let cmat = cmat.readonly().as_array().to_owned();
        Self {
            core: rswfn::base::Coefs::new(
                ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat,
            ),
        }
    }

    /// 由引用数据创建轨道系数对象
    ///
    /// **参数**
    /// - `ato_atms`: 每个基函数所属的原子索引
    /// - `ato_shls`: 每个基函数所属的壳层索引
    /// - `ato_syms`: 每个基函数的符号
    /// - `obt_engs`: 轨道能量列表
    /// - `obt_occs`: 轨道占据数列表
    /// - `cmat`: 轨道系数矩阵（NumPy 二维数组）
    #[staticmethod]
    pub fn new_ref(
        ato_atms: Vec<usize>,
        ato_shls: Vec<usize>,
        ato_syms: Vec<String>,
        obt_engs: Vec<f64>,
        obt_occs: Vec<f64>,
        cmat: Bound<'_, PyArray2<f64>>,
    ) -> Self {
        let cmat = cmat.readonly().as_array().to_owned();
        Self {
            core: rswfn::base::Coefs::new_ref(
                &ato_atms, &ato_shls, &ato_syms, &obt_engs, &obt_occs, &cmat,
            ),
        }
    }

    /// 判断轨道系数是否为空
    pub fn is_none(&self) -> bool {
        self.core.is_none()
    }

    /// 设置指定基函数形式的轨道系数矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式
    /// - `cmat`: 轨道系数矩阵（NumPy 二维数组）
    /// - `mole`: 分子对象
    pub fn set_cmat(&mut self, form: &str, cmat: Bound<'_, PyArray2<f64>>, mole: Mole) {
        let cmat = cmat.readonly().as_array().to_owned();
        self.core.set_cmat(form, &cmat, &mole.core)
    }

    /// 获取占据或虚轨道对应的电子数
    ///
    /// **参数**
    /// - `mole`: 分子对象
    ///
    /// **返回**
    /// 返回 [alpha 电子数, beta 电子数]
    pub fn obt_eles(&self, mole: Mole) -> [usize; 2] {
        self.core.obt_eles(&mole.core)
    }

    /// 判断两个轨道系数是否相同
    ///
    /// **参数**
    /// - `other`: 另一个轨道系数对象
    pub fn equal(&self, other: Coefs) -> bool {
        self.core.equal(&other.core)
    }

    /// 获取指定基函数形式的轨道系数矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    /// - `mole`: 分子对象
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_cmat(&self, py: Python, form: &str, mole: Mole) -> Py<PyArray2<f64>> {
        let cmat = self.core.get_cmat(form, &mole.core).to_owned();
        cmat.into_pyarray(py).unbind()
    }

    /// 获取指定基函数形式的密度矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    /// - `mole`: 分子对象
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_dmat(&self, py: Python, form: &str, mole: Mole) -> Py<PyArray2<f64>> {
        let dmat = self.core.get_dmat(form, &mole.core).to_owned();
        dmat.into_pyarray(py).unbind()
    }

    /// 获取所有轨道的占据数列表
    #[getter]
    pub fn obt_occs(&self) -> Vec<f64> {
        self.core.obt_occs.clone()
    }

    /// 获取所有轨道的能量本征值列表
    #[getter]
    pub fn obt_engs(&self) -> Vec<f64> {
        self.core.obt_engs.clone()
    }

    /// 获取占据或虚轨道的索引列表
    ///
    /// **参数**
    /// - `if_occ`: 若为 `true` 返回占据轨道索引，否则返回虚轨道索引
    pub fn obt_idxs(&self, if_occ: bool) -> Vec<usize> {
        self.core.obt_idxs(if_occ).clone()
    }

    /// 返回轨道系数信息的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
