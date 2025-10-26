use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::bondprop::direction::Calculator<'_> {
        rswfn::bondprop::direction::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn verts(&self, ai: u32, aj: u32) -> Vec<[f64; 3]> {
        let dirs = self.caler().verts(ai, aj);
        dirs.into_iter().map(|dir| [dir.x, dir.y, dir.z]).collect()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "direction")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
