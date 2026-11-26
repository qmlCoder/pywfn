use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct ObjWriter {
    core: rswfn::writer::ObjWriter,
}

#[pymethods]
impl ObjWriter {
    #[new]
    pub fn new(verts: Vec<[f64; 3]>, faces: Vec<[usize; 3]>) -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::ObjWriter::new(verts, faces),
        };
        Ok(writer)
    }

    #[getter(verts)]
    pub fn get_verts(&self) -> Vec<[f64; 3]> {
        self.core.verts.clone()
    }

    #[setter(verts)]
    pub fn set_verts(&mut self, verts: Vec<[f64; 3]>) {
        self.core.verts = verts;
    }

    #[getter(faces)]
    pub fn get_faces(&self) -> Vec<[usize; 3]> {
        self.core.faces.clone()
    }

    #[setter(faces)]
    pub fn set_faces(&mut self, faces: Vec<[usize; 3]>) {
        self.core.faces = faces;
    }

    pub fn save(&self, path: String) {
        self.core.save(&path);
    }
}
