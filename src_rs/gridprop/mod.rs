pub mod density;
pub mod dftgrid;
pub mod potential;
pub mod wfnfunc;

use pyo3::prelude::*;
use rswfn;
use rswfn::base::Vec3;

#[pyclass]
pub struct LineGrid {
    pub core: rswfn::gridprop::LineGrid,
}

#[pymethods]
impl LineGrid {
    #[new]
    pub fn new() -> Self {
        Self {
            core: rswfn::gridprop::LineGrid::new(),
        }
    }

    pub fn set(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64) {
        let p0 = Vec3::from_arr(p0);
        let p1 = Vec3::from_arr(p1);
        self.core.set(&p0, &p1, step);
    }

    pub fn get(&self) -> ([usize; 1], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct RectGrid {
    pub core: rswfn::gridprop::RectGrid,
}

#[pymethods]
impl RectGrid {
    #[new]
    pub fn new() -> Self {
        RectGrid {
            core: rswfn::gridprop::RectGrid::new(),
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
        let center = Vec3::from_arr(center);
        let normal = Vec3::from_arr(normal);
        let vector = Vec3::from_arr(vector);
        self.core.set_v1(&center, &normal, &vector, size, step);
    }

    pub fn set_v2(&mut self, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], step: [f64; 2]) {
        let p0 = Vec3::from_arr(p0);
        let p1 = Vec3::from_arr(p1);
        let p2 = Vec3::from_arr(p2);
        self.core.set_v2(&p0, &p1, &p2, step);
    }

    pub fn get(&self) -> ([usize; 2], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct CubeGrid {
    core: rswfn::gridprop::CubeGrid,
}

#[pymethods]
impl CubeGrid {
    #[new]
    pub fn new() -> Self {
        CubeGrid {
            core: rswfn::gridprop::CubeGrid::new(),
        }
    }

    pub fn set_v1(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64, bord: f64) {
        self.core.set_v1(&p0, &p1, step, bord);
    }

    pub fn get(&self) -> ([usize; 3], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
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

#[pymodule]
pub mod gridprop {
    #[pymodule_export]
    use super::CubeGrid;
    #[pymodule_export]
    use super::LineGrid;
    #[pymodule_export]
    use super::RectGrid;

    #[pymodule_export]
    use super::density::density;
    #[pymodule_export]
    use super::dftgrid::dftgrid;
    #[pymodule_export]
    use super::potential::potential;
    #[pymodule_export]
    use super::wfnfunc::wfnfunc;
}
