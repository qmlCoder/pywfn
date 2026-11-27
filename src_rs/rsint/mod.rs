#![allow(non_snake_case)]

use std::collections::HashMap;

use crate::base::{Mole, Stm};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct RsInt {
    pub mole: Mole,
}

impl RsInt {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::rsint::RsInt<'_> {
        rswfn::rsint::RsInt::new(&self.mole.core)
    }
}

#[pymethods]
impl RsInt {
    /// 创建一个新的 RsInt 积分计算实例
    ///
    /// ## 参数
    /// - `mole`: 包含分子、基组信息的 Mole 对象
    ///
    /// ## 返回
    /// 返回一个 RsInt 对象
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算重叠矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph"（球谐）或 "car"（笛卡尔）
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的重叠矩阵
    pub fn smat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let smat = self.caler().smat(&form).clone();
        smat.into_pyarray(py).unbind()
    }

    /// 计算动能矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的动能矩阵
    pub fn tmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let tmat = self.caler().tmat(&form).clone();
        tmat.into_pyarray(py).unbind()
    }

    /// 计算势能（核吸引）矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的势能矩阵
    pub fn vmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let vmat = self.caler().vmat(&form).clone();
        vmat.into_pyarray(py).unbind()
    }

    /// 同时计算重叠、动能和势能三个单电子矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `cals`: 指定哪些矩阵需要计算，[calc_s, calc_t, calc_v]，1 表示计算，0 表示跳过
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个包含三个二维 NumPy 数组的数组 [S, T, V]
    pub fn mat1e(&self, py: Python, cals: [usize; 3], form: String) -> [Py<PyArray2<f64>>; 3] {
        let [smat, tmat, vmat] = self.caler().mat1e(cals, &form);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }

    /// 在局域坐标系下同时计算三个单电子矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `stms`: 原子到壳层映射的字典（键为原子索引，值为 Stm 对象）
    /// - `cals`: 指定哪些矩阵需要计算，[calc_s, calc_t, calc_v]
    ///
    /// ## 返回
    /// 返回一个包含三个二维 NumPy 数组的数组 [S, T, V]
    pub fn mat1e_lcs(
        &self,
        py: Python,
        stms: HashMap<usize, Stm>,
        cals: [usize; 3],
    ) -> [Py<PyArray2<f64>>; 3] {
        let stms = stms.into_iter().map(|(idx, stm)| (idx, stm.core)).collect();
        let [smat, tmat, vmat] = self.caler().mat1e_lcs(&stms, cals);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsint")?; // 为当前的rs文件创建一个子模块
    m.add_class::<RsInt>()?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
