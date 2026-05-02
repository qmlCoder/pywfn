use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct CubWriter {
    core: rswfn::writer::CubWriter,
}

#[pymethods]
impl CubWriter {
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::CubWriter::new(),
        };
        Ok(writer)
    }

    pub fn set(
        &mut self,
        title: Option<String>,
        syms: Option<Vec<String>>,
        xyzs: Option<Vec<[f64; 3]>>,
        obts: Option<Vec<isize>>,
        pos0: Option<[f64; 3]>,
        size: Option<[usize; 3]>,
        step: Option<[f64; 3]>,
        vals: Option<Vec<Vec<f64>>>,
    ) {
        if let Some(title) = title {
            self.core.title = title
        }
        if let Some(syms) = syms {
            self.core.syms = syms
        }
        if let Some(xyzs) = xyzs {
            self.core.xyzs = xyzs;
        }
        if let Some(obts) = obts {
            self.core.obts = obts;
        }
        if let Some(pos0) = pos0 {
            self.core.pos0 = pos0;
        }
        if let Some(size) = size {
            self.core.size = size;
        }
        if let Some(step) = step {
            self.core.step = step;
        }
        if let Some(vals) = vals {
            self.core.vals = vals;
        }
    }

    pub fn read_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }

    pub fn save(&mut self, path: &str) {
        self.core.save(path);
    }
}
