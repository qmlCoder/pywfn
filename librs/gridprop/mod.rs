pub mod density;
pub mod dftgrid;
pub mod potential;
pub mod wfnfunc;

use nalgebra::Vector3;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct LineGrid {
    pub inner: rswfn::gridprop::LineGrid,
}

#[pymethods]
impl LineGrid {
    #[new]
    pub fn new() -> Self {
        Self {
            inner: rswfn::gridprop::LineGrid::new(),
        }
    }

    pub fn set(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64) {
        let p0 = Vector3::new(p0[0], p0[1], p0[2]);
        let p1 = Vector3::new(p1[0], p1[1], p1[2]);
        self.inner.set(p0, p1, step);
    }

    pub fn get(&self) -> ([u32; 1], Vec<[f64; 3]>) {
        let (shape, grids) = self.inner.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct RectGrid {
    pub inner: rswfn::gridprop::RectGrid,
}

#[pymethods]
impl RectGrid {
    #[new]
    pub fn new() -> Self {
        RectGrid {
            inner: rswfn::gridprop::RectGrid::new(),
        }
    }

    pub fn set_v1(
        &mut self,
        center: [f64; 3],
        normal: [f64; 3],
        vector: [f64; 3],
        size: [f64; 2],
        step: [f64; 2],
    ) {
        let center = Vector3::new(center[0], center[1], center[2]);
        let normal = Vector3::new(normal[0], normal[1], normal[2]);
        let vector = Vector3::new(vector[0], vector[1], vector[2]);
        self.inner.set_v1(&center, &normal, &vector, size, step);
    }

    pub fn set_v2(&mut self, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], step: [f64; 2]) {
        let p0 = Vector3::new(p0[0], p0[1], p0[2]);
        let p1 = Vector3::new(p1[0], p1[1], p1[2]);
        let p2 = Vector3::new(p2[0], p2[1], p2[2]);
        self.inner.set_v2(&p0, &p1, &p2, step);
    }

    pub fn get(&self) -> ([u32; 2], Vec<[f64; 3]>) {
        let (shape, grids) = self.inner.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct CubeGrid {
    inner: rswfn::gridprop::CubeGrid,
}

#[pymethods]
impl CubeGrid {
    #[new]
    pub fn new() -> Self {
        CubeGrid {
            inner: rswfn::gridprop::CubeGrid::new(),
        }
    }

    pub fn set_v1(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64, bord: f64) {
        self.inner.set_v1(&p0, &p1, step, bord);
    }

    pub fn get(&self) -> ([u32; 3], Vec<[f64; 3]>) {
        let (shape, grids) = self.inner.get();
        (shape, grids.clone())
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "gridprop")?;
    density::register_module(&m)?;
    dftgrid::register_module(&m)?;
    wfnfunc::register_module(&m)?;
    m.add_class::<LineGrid>()?;
    m.add_class::<RectGrid>()?;
    m.add_class::<CubeGrid>()?;
    parent_module.add_submodule(&m)
}
