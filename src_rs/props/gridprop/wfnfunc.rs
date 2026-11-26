use crate::base::Mole;
use crate::props::gridprop::vals2numpy;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3};
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    /// 创建底层的 `rswfn::props::gridprop::wfnfunc::Calculator` 实例。
    fn caler(&self) -> rswfn::props::gridprop::wfnfunc::Calculator<'_> {
        rswfn::props::gridprop::wfnfunc::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的波函数计算器。
    ///
    /// * `mole` - 分子对象
    pub fn new(mole: Mole) -> Self {
        let atms = mole.core.atoms().idxs().to_owned();
        Self { mole, atms }
    }

    pub fn ato_wfn(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> Vec<(Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>)> {
        let mut caler = self.caler();
        caler.atms = &self.atms;

        // 没有格点时直接返回空列表，避免后续 vall[0] 越界
        if grids.is_empty() {
            return Vec::new();
        }

        // [n_grids, n_aos]
        let vall: Vec<Vec<_>> = grids
            .into_par_iter()
            .map(|grid| caler.ato_wfn(&grid, level))
            .collect();

        let n_aos = vall[0].len();
        let n_grids = vall.len();

        // 转置为 [n_aos, n_grids]
        let mut result: Vec<Vec<_>> = (0..n_aos).map(|_| Vec::with_capacity(n_grids)).collect();
        for grid_vals in vall {
            for (iao, val) in grid_vals.into_iter().enumerate() {
                result[iao].push(val);
            }
        }

        result
            .into_iter()
            .map(|ao_vals| vals2numpy(py, &ao_vals, level))
            .collect()
    }

    pub fn obt_wfn(
        &self,
        py: Python,
        obt: usize,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let mut caler = self.caler();
        caler.atms = &self.atms;
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.obt_wfn(obt, &grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }
}

/// 计算高斯型轨道（Gaussian Type Orbital, GTO）在给定格点上的值。
///
/// * `xyzs` - 以原子为中心的笛卡尔坐标点列表
/// * `lmn` - 角动量指数 `[l, m, n]`，对应 `x^l * y^m * z^n`
/// * `alp` - 高斯指数参数
/// * `level` - 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）
///
/// 返回 `(波函数值, 梯度, Hessian)` 三元组。
#[pyfunction]
pub fn gtf(
    xyzs: Vec<[f64; 3]>,
    lmn: [usize; 3],
    alp: f64,
    level: usize,
) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>) {
    let vals: Vec<(f64, [f64; 3], [[f64; 3]; 3])> = xyzs
        .into_par_iter()
        .map(|xyz| rswfn::props::gridprop::wfnfunc::gtf(&xyz, &lmn, alp, level))
        .collect();
    let val0: Vec<f64> = vals.iter().map(|val| val.0).collect();
    let mut val1: Vec<[f64; 3]> = Vec::new();
    let mut val2: Vec<[[f64; 3]; 3]> = Vec::new();
    if level >= 1 {
        val1 = vals.iter().map(|val| val.1).collect();
    }
    if level >= 2 {
        val2 = vals.iter().map(|val| val.2).collect();
    }
    (val0, val1, val2)
}

/// 计算收缩高斯函数（Contracted Gaussian Function, CGF）在给定格点上的值。
///
/// CGF 是多个 GTO 的线性组合，用于表示原子轨道。
///
/// * `xyzs` - 以原子为中心的笛卡尔坐标点列表
/// * `lmn` - 角动量指数 `[l, m, n]`
/// * `alps` - 高斯指数参数列表
/// * `coes` - 收缩系数列表
/// * `level` - 计算级别（0=仅波函数值，1=值+梯度，2=值+梯度+Hessian）
///
/// 返回每个格点上的 `(波函数值, 梯度, Hessian)` 向量列表。
#[pyfunction]
pub fn cgf(
    xyzs: Vec<[f64; 3]>,
    lmn: [usize; 3],
    alps: Vec<f64>,
    coes: Vec<f64>,
    level: usize,
) -> Vec<(f64, [f64; 3], [[f64; 3]; 3])> {
    xyzs.into_par_iter()
        .map(|xyz| rswfn::props::gridprop::wfnfunc::cgf(&xyz, &lmn, &alps, &coes, level))
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "wfnfunc")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(gtf, &m)?)?;
    m.add_function(wrap_pyfunction!(cgf, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod wfnfunc {
    #[pymodule_export]
    use super::Calculator;
}
