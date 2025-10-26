use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::fragprop::energy::Calculator<'_> {
        rswfn::fragprop::energy::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }
    pub fn EIEBA(&self, fragA: Vec<u32>, fragB: Vec<u32>) -> (f64, f64, f64) {
        self.caler().EIEBA(&fragA, &fragB)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "energy")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
