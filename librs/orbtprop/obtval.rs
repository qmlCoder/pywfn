use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::orbtprop::obtval::Calculator<'_> {
        rswfn::orbtprop::obtval::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn pi_ele(&self) -> Vec<f64> {
        self.caler().pi_ele()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "obtval")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
