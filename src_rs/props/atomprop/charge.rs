use ndarray::{Array, Array1, Array2};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;
use rswfn::base::Vec3;

use crate::base::{Mole, Stm};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::atomprop::charge::Calculator<'_> {
        let mut caler = rswfn::props::atomprop::charge::Calculator::new(&self.mole.core);
        caler.atms = &self.atms;
        caler
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的电荷计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `atms`: 可选的需要计算的原子索引列表，默认为所有原子
    pub fn new(mole: Mole, atms: Option<Vec<usize>>) -> Self {
        let atms = match atms {
            Some(atms) => atms,
            None => mole.core.atoms().idxs().clone(),
        };
        Self { mole, atms }
    }

    /// 计算自旋密度
    ///
    /// 返回每个原子的自旋密度分布矩阵，包含 alpha 和 beta 电子密度差。
    ///
    /// ## 参数
    /// - `py`: Python 解释器引用
    /// - `ctype`: 电荷类型
    ///
    /// ## 返回
    /// N×3 的 NumPy 二维数组，每行对应一个原子的自旋密度分量
    pub fn spin(&self, py: Python, ctype: &str) -> Py<PyArray2<f64>> {
        let vals = self.caler().spin(ctype);
        let arr = Array2::from_shape_vec((vals.len(), 3), vals.concat()).unwrap();
        arr.into_pyarray(py).unbind()
    }

    /// 计算 Mulliken 布居分析电荷
    ///
    /// 基于 Mulliken 方法计算每个原子的净电荷。
    ///
    /// ## 参数
    /// - `py`: Python 解释器引用
    ///
    /// ## 返回
    /// 一维 NumPy 数组，每个元素对应一个原子的 Mulliken 电荷
    pub fn mulliken(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().mulliken();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    /// 计算 Löwdin 布居分析电荷
    ///
    /// 基于 Löwdin 方法计算每个原子的净电荷。
    ///
    /// ## 参数
    /// - `py`: Python 解释器引用
    ///
    /// ## 返回
    /// 一维 NumPy 数组，每个元素对应一个原子的 Löwdin 电荷
    pub fn lowdin(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().lowdin();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    /// 计算 Hirshfeld 布居分析电荷
    ///
    /// 基于 Hirshfeld 方法计算每个原子的净电荷。
    ///
    /// ## 参数
    /// - `py`: Python 解释器引用
    ///
    /// ## 返回
    /// 一维 NumPy 数组，每个元素对应一个原子的 Hirshfeld 电荷
    pub fn hirshfeld(&self, py: Python) -> Py<PyArray1<f64>> {
        let vals = self.caler().hirshfeld();
        let arr = Array1::from_vec(vals);
        arr.into_pyarray(py).unbind()
    }

    /// 计算定域分子轨道上的电荷分布
    ///
    /// 基于指定的原子轨道方向，计算定域分子轨道中的电荷分布信息。
    ///
    /// ## 参数
    /// - `dirs`: 原子索引到方向向量的映射，指定每个原子的轨道方向
    /// - `keep_other_atm`: 是否保留其他原子的贡献
    /// - `keep_other_sym`: 是否保留其他对称性的贡献
    /// - `ctype`: 电荷类型
    ///
    /// ## 返回
    /// 每个定域轨道的电荷值列表
    pub fn pocv(
        &self,
        dirs: HashMap<usize, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let dirs: HashMap<usize, Vec3> = dirs
            .into_iter()
            .map(|(atm, dir)| (atm, Vec3::from_arr(dir)))
            .collect();
        let vals = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        Ok(vals)
    }

    /// 计算分子轨道上的电荷分布
    ///
    /// 基于指定的原子轨道旋转矩阵，计算分子轨道中的电荷分布。
    ///
    /// ## 参数
    /// - `stms`: 原子索引到旋转矩阵（Stm）的映射
    /// - `atos`: 需要计算的原子轨道索引列表
    /// - `ctype`: 电荷类型
    ///
    /// ## 返回
    /// 每个分子轨道的电荷值列表
    pub fn mocv(
        &self,
        stms: HashMap<usize, Stm>,
        atos: Vec<usize>,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let stms = stms
            .into_iter()
            .map(|(atm, stm)| (atm, stm.core.clone()))
            .collect();
        let vals = self.caler().mocv(&stms, &atos, ctype);
        Ok(vals)
    }

    /// 计算 π 定域分子轨道上的电荷分布
    ///
    /// 提取 π 电子体系的定域分子轨道电荷分布，包含 π 轨道方向信息。
    ///
    /// ## 参数
    /// - `ctype`: 电荷类型
    ///
    /// ## 返回
    /// - `HashMap<usize, [f64; 3]>`: 原子索引到 π 轨道方向的映射
    /// - `Vec<f64>`: 每个 π 定域轨道的电荷值列表
    pub fn pi_pocv(&self, ctype: &str) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<f64>)> {
        let (dirs, eles) = self
            .caler()
            .pi_pocv(ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, dir.into_arr()))
            .collect();
        Ok((dirs, eles))
    }

    pub fn pi_mocv(&self, ctype: &str) -> PyResult<(HashMap<usize, Stm>, Vec<usize>, Vec<f64>)> {
        let (stms, atos, eles) = self
            .caler()
            .pi_mocv(ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let stms = stms
            .iter()
            .map(|(atm, stm)| (*atm, Stm { core: stm.clone() }))
            .collect();
        Ok((stms, atos, eles))
    }

    pub fn pi_mocv_v2(&self, ctype: &str) -> PyResult<(HashMap<usize, Stm>, Vec<usize>, Vec<f64>)> {
        let (stms, atos, eles) = self
            .caler()
            .pi_mocv_v2(ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let stms = stms
            .iter()
            .map(|(atm, stm)| (*atm, Stm { core: stm.clone() }))
            .collect();
        Ok((stms, atos, eles))
    }
}

/// 计算 Mulliken 布居分析（独立函数）
///
/// 直接使用密度矩阵和重叠矩阵计算 Mulliken 布居分析。
///
/// ## 参数
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `smat`: 重叠矩阵（NumPy 二维数组）
/// - `uplo`: 上三角/下三角的原子轨道索引对列表
///
/// ## 返回
/// 每个原子的 Mulliken 电荷列表
#[pyfunction]
pub fn mulliken(
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Vec<f64> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::props::atomprop::charge::mulliken(&dmat, &smat, &uplo)
}

/// 计算 Löwdin 布居分析（独立函数）
///
/// 直接使用密度矩阵和重叠矩阵计算 Löwdin 布居分析。
///
/// ## 参数
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `smat`: 重叠矩阵（NumPy 二维数组）
/// - `uplo`: 上三角/下三角的原子轨道索引对列表
///
/// ## 返回
/// 每个原子的 Löwdin 电荷列表
#[pyfunction]
pub fn lowdin(
    dmat: Bound<'_, PyArray2<f64>>,
    smat: Bound<'_, PyArray2<f64>>,
    uplo: Vec<[usize; 2]>,
) -> Vec<f64> {
    let dmat = dmat.readonly().as_array().to_owned();
    let smat = smat.readonly().as_array().to_owned();
    rswfn::props::atomprop::charge::lowdin(&dmat, &smat, &uplo)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "charge")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(mulliken, &m)?)?;
    m.add_function(wrap_pyfunction!(lowdin, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod charge {
    #[pymodule_export]
    use super::Calculator;
}
