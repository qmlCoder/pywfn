use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct FchWriter {
    core: rswfn::writer::FchWriter,
}

#[pymethods]
impl FchWriter {
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::FchWriter::new(),
        };
        Ok(writer)
    }

    #[staticmethod]
    pub fn from_mole(mole: Mole) -> Self {
        Self {
            core: rswfn::writer::FchWriter::from_mole(&mole.core),
        }
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

    pub fn read_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }
}
