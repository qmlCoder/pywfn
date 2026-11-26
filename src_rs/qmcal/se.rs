use crate::base::{Coefs, Mole};
use pyo3::prelude::*;
use rswfn;

/// RSE（限制性开壳层自洽场）计算器
///
/// 提供开壳层体系的扩展 Hückel 近似计算功能。
#[pyclass]
pub struct RSE {
    mole: Mole,
}

impl RSE {
    fn caler(&self) -> rswfn::qmcal::se::RSE<'_> {
        rswfn::qmcal::se::RSE::new(&self.mole.core)
    }
}

#[pymethods]
impl RSE {
    /// 创建一个新的 RSE 计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> PyResult<Self> {
        let rse = Self { mole };
        Ok(rse)
    }

    /// 执行扩展 Hückel 计算。
    ///
    /// 对分子进行扩展 Hückel 近似计算，返回总能量和分子轨道系数。
    ///
    /// # 返回
    /// * `(f64, Coefs)` - 元组，包含计算得到的总能量和分子轨道系数。
    pub fn run(&self) -> (f64, Coefs) {
        let (eng, core) = self.caler().calc_rse(false);
        (eng, Coefs { core })
    }
}

/// 在父模块中注册 `se` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "se")?;
    m.add_class::<RSE>()?;
    parent_module.add_submodule(&m)
}
