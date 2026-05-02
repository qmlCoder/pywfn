use crate::base::Mole;
use ndarray::Array4;
use numpy::PyReadonlyArray4;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct RHF {
    mole: Mole,
}

impl RHF {
    fn caler(&self) -> rswfn::qmcalc::hf::RHF<'_> {
        rswfn::qmcalc::hf::RHF::new(&self.mole.core)
    }
}

#[pymethods]
impl RHF {
    #[new]
    pub fn new(mole: Mole) -> PyResult<Self> {
        let rhf = Self { mole };
        Ok(rhf)
    }

    pub fn run(&self, emat: Option<PyReadonlyArray4<f64>>) {
        let emat = match emat {
            Some(emat) => {
                let emat = emat.as_array().to_owned();
                let emat = rswfn::integ::Mat2e::from_mat4(emat);
                Some(emat)
            }
            None => None,
        };
        self.caler().run(emat);
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "hf")?;
    m.add_class::<RHF>()?;
    parent_module.add_submodule(&m)
}
