use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::moleprop::aromacity::Calculator<'_> {
        rswfn::moleprop::aromacity::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 使用pi键级标准差计算芳香性
    pub fn PISD(&self, rings: Option<Vec<Vec<u32>>>) -> Vec<f64> {
        self.caler().PISD(rings.as_ref())
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "aromacity")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
