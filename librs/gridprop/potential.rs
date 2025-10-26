use crate::base::Mole;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::gridprop::potential::Calculator<'_> {
        rswfn::gridprop::potential::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn nuc_potential(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.nuc_potential(&grid))
            .collect()
    }

    pub fn ele_potential(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.ele_potential(&grid))
            .collect()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "potential")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
