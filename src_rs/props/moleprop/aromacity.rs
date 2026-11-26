use crate::base::Mole;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// 芳香性计算器
///
/// 基于分子轨道信息计算分子的芳香性指数。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::moleprop::aromacity::Calculator<'_> {
        rswfn::props::moleprop::aromacity::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的芳香性计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 使用 π 键级的标准差计算芳香性指数。
    ///
    /// 对于每个指定的环，计算环上原子间 π 键级的标准差，
    /// 标准差越小表明该环的芳香性越强。
    ///
    /// # 参数
    /// * `rings` - 可选，指定要计算的环列表，每个环由原子索引列表表示。
    ///             如果为 `None`，则自动检测所有环。
    ///
    /// # 返回
    /// * `Vec<f64>` - 每个环对应的芳香性指数值列表。
    pub fn PISD(&self, rings: Option<Vec<Vec<usize>>>) -> PyResult<Vec<f64>> {
        self.caler()
            .PISD(rings.as_ref())
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))
    }
}

/// 在父模块中注册 `aromacity` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "aromacity")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod aromacity {
    #[pymodule_export]
    use super::Calculator;
}
