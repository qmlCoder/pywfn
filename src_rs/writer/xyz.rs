use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct XyzWriter {
    core: rswfn::writer::XyzWriter,
}

#[pymethods]
impl XyzWriter {
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::XyzWriter::new(),
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

    #[getter(title)]
    pub fn get_title(&self) -> String {
        self.core.title.clone()
    }

    #[setter(title)]
    pub fn set_title(&mut self, title: String) {
        self.core.title = title;
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
