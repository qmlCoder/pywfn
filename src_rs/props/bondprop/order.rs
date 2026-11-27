use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use crate::base::{Mole, Stm};
// use crate::props::orbtprop::obtmat::Mocv;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::bondprop::order::Calculator<'_> {
        rswfn::props::bondprop::order::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的键级计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mayer(&self, py: Python) -> Py<PyArray2<f64>> {
        let omat = self.caler().mayer();
        let omat = omat.into_pyarray(py).unbind();
        omat
    }

    pub fn wiberg(&self, py: Python) -> Py<PyArray2<f64>> {
        let omat = self.caler().wiberg();
        let omat = omat.into_pyarray(py).unbind();
        omat
    }

    pub fn order_obt(&self, py: Python, otype: &str, iobt: usize) -> Py<PyArray2<f64>> {
        let omat = self.caler().order_obt(otype, iobt);
        let omat = omat.into_pyarray(py).unbind();
        omat
    }

    pub fn pocv(
        &self,
        py: Python,
        dirs: HashMap<usize, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        otype: &str,
    ) -> PyResult<Py<PyArray2<f64>>> {
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, Vec3::from_arr(dir)))
            .collect();
        let omat = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, otype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        Ok(omat.into_pyarray(py).unbind())
    }

    pub fn mocv(
        &self,
        py: Python,
        stms: HashMap<usize, Stm>, // 每个原子的保留信息
        atos: Vec<usize>,
        otype: &str,
    ) -> Py<PyArray2<f64>> {
        let stms = stms.into_iter().map(|(key, stm)| (key, stm.core)).collect();
        let omat = self.caler().mocv(&stms, &atos, otype);
        omat.into_pyarray(py).unbind()
    }

    pub fn pi_pocv(
        &self,
        py: Python,
        otype: &str,
    ) -> PyResult<(HashMap<usize, [f64; 3]>, Py<PyArray2<f64>>)> {
        let (dirs, omat) = self
            .caler()
            .pi_pocv(otype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let dirs = dirs
            .into_iter()
            .map(|(key, dir)| (key, dir.into_arr()))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        Ok((dirs, omat))
    }

    pub fn pi_mocv(
        &self,
        py: Python,
        otype: &str,
    ) -> (HashMap<usize, Stm>, Vec<usize>, Py<PyArray2<f64>>) {
        let (stms, atos, omat) = self.caler().pi_mocv(otype);
        let stms = stms
            .into_iter()
            .map(|(key, stm)| (key, Stm { core: stm }))
            .collect();
        let omat = omat.into_pyarray(py).unbind();
        (stms, atos, omat)
    }

    pub fn decom(&self, ai: usize, aj: usize, otype: String) -> Vec<f64> {
        self.caler().decom(ai, aj, &otype)
    }

    /// 计算指定方向的键级
    ///
    /// 在给定原子上沿特定方向计算与相邻原子间的键级。
    ///
    /// ## 参数
    /// - `atm`: 原子索引
    /// - `dir`: 方向向量 [x, y, z]
    /// - `nebs`: 相邻原子索引列表
    /// - `otype`: 轨道类型（如 "sigma", "pi"）
    ///
    /// ## 返回
    /// 与每个相邻原子在该方向上的键级值列表
    pub fn bound(
        &self,
        atm: usize,
        dir: [f64; 3],
        nebs: Vec<usize>,
        otype: &str,
        mtype: &str,
    ) -> PyResult<Vec<f64>> {
        let dir = Vec3::from_arr(dir);
        // rswfn 已将该方法拆分为 pocv/mocv，默认沿用旧实现的 pocv
        self.caler()
            .bound(atm, &dir, &nebs, otype, mtype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))
    }
    /// 按指定键级类型计算原子轨道键级矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `otype`: 键级类型，如 "mayer"、"wiberg"
    ///
    /// ## 返回
    /// NumPy 二维数组
    pub fn order(&self, py: Python, otype: &str) -> Py<PyArray2<f64>> {
        self.caler().order(otype).into_pyarray(py).unbind()
    }
}

/// 按指定键级类型计算原子轨道键级矩阵（独立函数）
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `smat`: 重叠矩阵（NumPy 二维数组）
/// - `otype`: 键级类型
///
/// ## 返回
/// NumPy 二维数组
#[pyfunction]
#[pyo3(name = "order")]
pub fn order_mat(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
    otype: &str,
) -> Py<PyArray2<f64>> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::props::bondprop::order::order(&dmat, &smat, otype)
        .into_pyarray(py)
        .unbind()
}

/// 计算 Mayer 键级（独立函数）
///
/// 直接使用密度矩阵和重叠矩阵计算 Mayer 键级矩阵。
///
/// ## 参数
/// - `py`: Python 解释器引用
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `smat`: 重叠矩阵（NumPy 二维数组）
///
/// ## 返回
/// N×N 的 NumPy 二维数组，表示 Mayer 键级矩阵
#[pyfunction]
pub fn mayer(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
) -> Py<PyArray2<f64>> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::props::bondprop::order::mayer(&dmat, &smat)
        .into_pyarray(py)
        .unbind()
}

/// 计算 Wiberg 键级（独立函数）
///
/// 直接使用密度矩阵和重叠矩阵计算 Wiberg 键级矩阵。
///
/// ## 参数
/// - `py`: Python 解释器引用
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `smat`: 重叠矩阵（NumPy 二维数组）
///
/// ## 返回
/// N×N 的 NumPy 二维数组，表示 Wiberg 键级矩阵
#[pyfunction]
pub fn wiberg(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
) -> Py<PyArray2<f64>> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::props::bondprop::order::wiberg(&dmat, &smat)
        .into_pyarray(py)
        .unbind()
}

/// 将原子轨道键级矩阵转换为原子间键级矩阵
///
/// 将原子轨道级别的键级矩阵按原子分组，得到原子间的总键级矩阵。
///
/// ## 参数
/// - `py`: Python 解释器引用
/// - `omat`: 原子轨道级别的键级矩阵（NumPy 二维数组）
/// - `uplo`: 上三角/下三角的原子轨道索引对列表
///
/// ## 返回
/// N×N 的 NumPy 二维数组，表示原子间的键级矩阵
#[pyfunction]
pub fn ato2atm(
    py: Python,
    omat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Py<PyArray2<f64>> {
    let omat = omat.readonly().as_array().to_owned();
    rswfn::props::bondprop::order::ato2atm(&omat, &uplo)
        .into_pyarray(py)
        .unbind()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "order")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(mayer, &m)?)?;
    m.add_function(wrap_pyfunction!(wiberg, &m)?)?;
    m.add_function(wrap_pyfunction!(order_mat, &m)?)?;
    m.add_function(wrap_pyfunction!(ato2atm, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod order {
    #[pymodule_export]
    use super::Calculator;
}
