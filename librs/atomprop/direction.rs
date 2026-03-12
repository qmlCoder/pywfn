use pyo3::prelude::*;

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

    pub fn normal_vector(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().normal_vector(atm);
        match dir {
            Some(dir) => Some([dir.x, dir.y, dir.z]),
            None => None,
        }
    }

    pub fn get_normal(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().get_normal(atm);
        match dir {
            Some(dir) => Some([dir.x, dir.y, dir.z]),
            None => None,
        }
    }

    pub fn LCS(&mut self, atm: usize, neb: Option<usize>) -> Option<[[f64; 3]; 3]> {
        let stm = self.caler().LCS(atm, neb);
        match stm {
            Some(stm) => {
                let [xx, xy, xz] = stm.ex;
                let [yx, yy, yz] = stm.ey;
                let [zx, zy, zz] = stm.ez;
                Some([[xx, xy, xz], [yx, yy, yz], [zx, zy, zz]])
            }
            None => None,
        }
    }

    pub fn reactions(&self, atm: usize) -> Vec<[f64; 3]> {
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
