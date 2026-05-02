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

    pub fn set_syms(&mut self, syms: Vec<String>) {
        self.core.syms = syms;
    }

    pub fn set_xyzs(&mut self, xyzs: Vec<[f64; 3]>) {
        self.core.xyzs = xyzs;
    }

    pub fn set_chk(&mut self, chk: String) {
        self.core.chk = chk;
    }

    pub fn set_job(&mut self, job: String) {
        self.core.job = job;
    }

    pub fn set_charge(&mut self, charge: isize) {
        self.core.charge = charge;
    }

    pub fn set_multip(&mut self, multip: usize) {
        self.core.multip = multip;
    }

    pub fn from_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }

    pub fn build(&self) -> String {
        self.core.build()
    }

    pub fn save(&self, path: String) {
        self.core.save(&path);
    }
}
