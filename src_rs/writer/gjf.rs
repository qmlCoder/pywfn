use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct GjfWriter {
    core: rswfn::writer::GjfWriter,
}

#[pymethods]
impl GjfWriter {
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::GjfWriter::new(),
        };
        Ok(writer)
    }

    #[getter(syms)]
    pub fn get_syms(&self) -> Vec<String> {
        self.core.syms.clone()
    }

    #[setter(syms)]
    pub fn set_syms(&mut self, syms: Vec<String>) {
        self.core.syms = syms;
    }

    #[getter(xyzs)]
    pub fn get_xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs.clone()
    }

    #[setter(xyzs)]
    pub fn set_xyzs(&mut self, xyzs: Vec<[f64; 3]>) {
        self.core.xyzs = xyzs;
    }

    #[getter(chk)]
    pub fn get_chk(&self) -> String {
        self.core.chk.clone()
    }

    #[setter(chk)]
    pub fn set_chk(&mut self, chk: String) {
        self.core.chk = chk;
    }

    #[getter(job)]
    pub fn get_job(&self) -> String {
        self.core.job.clone()
    }

    #[setter(job)]
    pub fn set_job(&mut self, job: String) {
        self.core.job = job;
    }

    #[getter(charge)]
    pub fn get_charge(&self) -> isize {
        self.core.charge
    }

    #[setter(charge)]
    pub fn set_charge(&mut self, charge: isize) {
        self.core.charge = charge;
    }

    #[getter(multip)]
    pub fn get_multip(&self) -> usize {
        self.core.multip
    }

    #[setter(multip)]
    pub fn set_multip(&mut self, multip: usize) {
        self.core.multip = multip;
    }

    pub fn build(&self) -> String {
        self.core.build()
    }

    pub fn save(&self, path: String) {
        self.core.save(&path);
    }
}
