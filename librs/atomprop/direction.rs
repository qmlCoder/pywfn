use pyo3::prelude::*;
use std::collections::HashMap;

use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::atomprop::direction::Calculator<'_> {
        rswfn::atomprop::direction::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn normal_vector(&self, atm: u32) -> Option<[f64; 3]> {
        let dir = self.caler().normal_vector(atm);
        match dir {
            Some(dir) => Some([dir.x, dir.y, dir.z]),
            None => None,
        }
    }

    pub fn get_normal(&self, atm: u32) -> Option<[f64; 3]> {
        let dir = self.caler().get_normal(atm);
        match dir {
            Some(dir) => Some([dir.x, dir.y, dir.z]),
            None => None,
        }
    }

    pub fn local_coord_system(&mut self, atm: u32, neb: Option<u32>) -> Option<[[f64; 3]; 3]> {
        let lcs = self.caler().local_coord_system(atm, neb);
        match lcs {
            Some(lcs) => {
                let xx = lcs[(0, 0)];
                let xy = lcs[(1, 0)];
                let xz = lcs[(2, 0)];

                let yx = lcs[(0, 1)];
                let yy = lcs[(1, 1)];
                let yz = lcs[(2, 1)];

                let zx = lcs[(0, 2)];
                let zy = lcs[(1, 2)];
                let zz = lcs[(2, 2)];
                Some([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])
            }
            None => None,
        }
    }

    pub fn reactions(&self, atm: u32) -> Vec<[f64; 3]> {
        let dirs = self.caler().reactions(atm);
        let dirs: Vec<[f64; 3]> = dirs.iter().map(|dir| [dir.x, dir.y, dir.z]).collect();
        dirs
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "direction")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
