use crate::base::Mole;
use nalgebra::Vector3;
use pyo3::prelude::*;
use rswfn;
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
            .map(|(key, dir)| (*key, [dir.x, dir.y, dir.z]))
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
        let dir = Vector3::new(dir[0], dir[1], dir[2]);
        let vals = self
            .caler()
            .fukui_dir(atm, &dir, &mole_n.core, &mole_p.core, ctype);
        vals
    }

    pub fn vector(&self, atm: usize, dir: [f64; 3], otype: &str) -> (Vec<usize>, Vec<f64>, f64) {
        let dir = Vector3::new(dir[0], dir[1], dir[2]);
        self.caler().vector(atm, &dir, otype)
    }

    pub fn free_valence(&self) -> Vec<f64> {
        self.caler().free_valence()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "activity")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
