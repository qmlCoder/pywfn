use pyo3::prelude::*;

use rswfn;

use numpy::{IntoPyArray, PyArray2};

use crate::base::Mole;

/// 布居分析计算器
///
/// 提供分子轨道布居分析相关的计算功能。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::moleprop::population::Calculator<'_> {
        rswfn::props::moleprop::population::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的布居分析计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }
}

/// 在父模块中注册 `population` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "population")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod population {
    #[pymodule_export]
    use super::Calculator;
}
