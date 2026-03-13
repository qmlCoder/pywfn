use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

use crate::base::Mole;
use crate::orbtprop::obtmat::Mocv;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::charge::Calculator<'_> {
        let mut caler = rswfn::atomprop::charge::Calculator::new(&self.mole.core);
        caler.atms = &self.atms;
        caler
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole, atms: Option<Vec<usize>>) -> Self {
        let atms = match atms {
            Some(atms) => atms,
            None => mole.core.atoms().idxs().clone(),
        };
        Self { mole, atms }
    }

    pub fn spin(&self, ctype: &str) -> Vec<[f64; 3]> {
        self.caler().spin(ctype)
    }

    pub fn mulliken(&self) -> Vec<f64> {
        self.caler().mulliken()
    }

    pub fn lowdin(&self) -> Vec<f64> {
        self.caler().lowdin()
    }

    pub fn hirshfeld(&self) -> Vec<f64> {
        self.caler().hirshfeld()
    }

    pub fn pocv(
        &self,
        dirs: HashMap<usize, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let dirs: HashMap<usize, Vector3<f64>> = dirs
            .into_iter()
            .map(|(atm, dir)| {
                let [x, y, z] = dir;
                (atm, Vector3::new(x, y, z))
            })
            .collect();
        let vals = self
            .caler()
            .pocv(&dirs, keep_other_atm, keep_other_sym, ctype);
        Ok(vals)
    }

    pub fn mocv(&self, mocvs: HashMap<usize, Mocv>, ctype: &str) -> PyResult<Vec<f64>> {
        let mocvs = mocvs
            .into_iter()
            .map(|(atm, mocv)| (atm, mocv.core.clone()))
            .collect();
        let vals = self.caler().mocv(&mocvs, ctype);
        Ok(vals)
    }

    pub fn pi_pocv(&self, ctype: &str) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<f64>)> {
        let (dirs, eles) = self.caler().pi_pocv(ctype);
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, [dir.x, dir.y, dir.z]))
            .collect();
        Ok((dirs, eles))
    }

    pub fn pi_mocv(&self, ctype: &str) -> PyResult<(HashMap<usize, Mocv>, Vec<f64>)> {
        let (mocvs, eles) = self.caler().pi_mocv(ctype);
        let mocvs = mocvs
            .iter()
            .map(|(atm, mocv)| (*atm, Mocv { core: mocv.clone() }))
            .collect();
        Ok((mocvs, eles))
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "charge")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
