use crate::base::{Coefs, Mole};
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct RSE {
    mole: Mole,
}

impl RSE {
    fn caler(&self) -> rswfn::qmcalc::se::RSE<'_> {
        rswfn::qmcalc::se::RSE::new(&self.mole.core)
    }
}

#[pymethods]
impl RSE {
    #[new]
    pub fn new(mole: Mole) -> PyResult<Self> {
        let rse = Self { mole };
        Ok(rse)
    }

    pub fn run(&self) -> (f64, Coefs) {
        let (eng, core) = self.caler().calc_rse(false);
        (eng, Coefs { core })
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "se")?;
    m.add_class::<RSE>()?;
    parent_module.add_submodule(&m)
}
