use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::moleprop::newmole::Calculator<'_> {
        rswfn::moleprop::newmole::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "newmole")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
