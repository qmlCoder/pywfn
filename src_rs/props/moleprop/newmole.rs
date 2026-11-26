use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

/// 新建分子计算器
///
/// 提供新建分子相关属性的计算功能。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::moleprop::variants::Calculator<'_> {
        rswfn::props::moleprop::variants::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的新建分子计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }
}

/// 在父模块中注册 `newmole` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "newmole")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod newmole {
    #[pymodule_export]
    use super::Calculator;
}
