use nalgebra::Vector3;
use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

use crate::base::Mole;
use crate::orbtprop::obtmat::Deco;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::charge::Calculator<'_> {
        rswfn::atomprop::charge::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn spin(&self) -> Vec<f64> {
        self.caler().spin()
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
        dirs: HashMap<u32, [f64; 3]>,
        keep_other_atm: bool,
        keep_other_sym: bool,
        ctype: &str,
    ) -> PyResult<Vec<f64>> {
        let dirs: HashMap<u32, Vector3<f64>> = dirs
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

    pub fn deco(&self, decos: HashMap<u32, Deco>, ctype: &str) -> PyResult<Vec<f64>> {
        let decos = decos
            .into_iter()
            .map(|(atm, deco)| (atm, deco.inner.clone()))
            .collect();
        let vals = self.caler().deco(&decos, ctype);
        Ok(vals)
    }

    pub fn pi_pocv(
        &self,
        ctype: &str,
        atms: Vec<u32>,
    ) -> PyResult<(HashMap<u32, [f64; 3]>, Vec<f64>)> {
        let (dirs, eles) = self.caler().pi_pocv(ctype, &atms);
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, [dir.x, dir.y, dir.z]))
            .collect();
        Ok((dirs, eles))
    }

    pub fn pi_deco(&self, ctype: &str) -> PyResult<(HashMap<u32, Deco>, Vec<f64>)> {
        let (decos, eles) = self.caler().pi_deco(ctype);
        let decos = decos
            .iter()
            .map(|(atm, deco)| {
                (
                    *atm,
                    Deco {
                        inner: deco.clone(),
                    },
                )
            })
            .collect();
        Ok((decos, eles))
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "charge")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
