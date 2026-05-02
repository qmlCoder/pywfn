use crate::base::Mole;
use pyo3::prelude::*;
use rswfn::{self, base::Vec3};
use std::collections::HashMap;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::activity::Calculator<'_> {
        rswfn::atomprop::activity::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn fukui(&self, mole_n: Mole, mole_p: Mole, ctype: &str) -> PyResult<Vec<[f64; 7]>> {
        let vals = self.caler().fukui(&mole_n.core, &mole_p.core, ctype);
        Ok(vals)
    }

    pub fn fukui_pi(
        &self,
        mole_n: Mole,
        mole_p: Mole,
        ctype: &str,
    ) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<[f64; 7]>)> {
        let (dirs, vals) = self.caler().fukui_pi(&mole_n.core, &mole_p.core, ctype);
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, dir.into_arr()))
            .collect();
        Ok((dirs, vals))
    }

    pub fn fukui_dir(
        &self,
        atm: usize,
        dir: [f64; 3],
        mole_n: Mole,
        mole_p: Mole,
        ctype: &str,
    ) -> [f64; 7] {
        let dir = Vec3::from_arr(dir);
        let vals = self
            .caler()
            .fukui_dir(atm, &dir, &mole_n.core, &mole_p.core, ctype);
        vals
    }

    pub fn freev(&self, otype: &str) -> Vec<[f64; 3]> {
        self.caler().freev(otype)
    }

    pub fn freev_pi(&self, otype: &str) -> (HashMap<usize, [f64; 3]>, Vec<[f64; 3]>) {
        let (dirs, vals) = self.caler().freev_pi(otype);
        let dirs = dirs
            .iter()
            .map(|(atm, dir)| (*atm, dir.into_arr()))
            .collect();
        (dirs, vals)
    }

    pub fn freev_dir(
        &self,
        atm: usize,
        dir: [f64; 3],
        nebs: Vec<usize>,
        otype: &str,
    ) -> (f64, Vec<f64>, f64) {
        let dir = Vec3::from_arr(dir);
        self.caler().freev_dir(atm, &dir, &nebs, otype)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "activity")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod activity {
    #[pymodule_export]
    use super::Calculator;
}
