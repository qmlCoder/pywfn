use pyo3::prelude::*;

use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use rswfn;

use crate::{base::Mole, utils};

#[derive(Clone)]
#[pyclass]
pub struct Coefs {
    pub core: rswfn::base::Coefs,
}

#[pymethods]
impl Coefs {
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
