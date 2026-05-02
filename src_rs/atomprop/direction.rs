use pyo3::prelude::*;

use rswfn;

use crate::base::{Mole, Stm};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::direction::Calculator<'_> {
        rswfn::atomprop::direction::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn normal_vector(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().normal_vector(atm);
        match dir {
            Some(dir) => Some(dir.into_arr()),
            None => None,
        }
    }

    pub fn get_normal(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().get_normal(atm);
        match dir {
            Some(dir) => Some(dir.into_arr()),
            None => None,
        }
    }

    pub fn LCS(&mut self, atm: usize, neb: Option<usize>) -> Option<Stm> {
        let stm = self.caler().LCS(atm, neb);
        match stm {
            Some(stm) => Some(Stm { core: stm }),
            None => None,
        }
    }

    pub fn reactions(&self, atm: usize) -> Vec<[f64; 3]> {
        let dirs = self.caler().reactions(atm);
        let dirs: Vec<[f64; 3]> = dirs.iter().map(|dir| dir.into_arr()).collect();
        dirs
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "direction")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod direction {
    #[pymodule_export]
    use super::Calculator;
}
