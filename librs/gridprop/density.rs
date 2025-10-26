use crate::base::Mole;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::gridprop::density::Calculator<'_> {
        rswfn::gridprop::density::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mol_rho(
        &self,
        grids: Vec<[f64; 3]>,
        level: u32,
    ) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>) {
        let caler = self.caler();
        let vals: Vec<(f64, [f64; 3], [[f64; 3]; 3])> = grids
            .into_par_iter()
            .map(|grid| caler.mol_rho_dm(&grid, level, None))
            .collect();
        let rho0: Vec<f64> = vals.iter().map(|val| val.0).collect();
        let rho1 = if level > 0 {
            vals.iter().map(|val| val.1).collect()
        } else {
            Vec::new()
        };
        let rho2 = if level > 0 {
            vals.iter().map(|val| val.2).collect()
        } else {
            Vec::new()
        };
        (rho0, rho1, rho2)
    }

    pub fn pi_rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.pi_rho(&grid))
            .collect()
    }

    pub fn ato_rho(
        &self,
        grids: Vec<[f64; 3]>,
        level: u32,
    ) -> Vec<(Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>)> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.ato_rho(&grid, level))
            .collect()
    }

    pub fn pro_rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.pro_rho(&grid))
            .collect()
    }

    pub fn obt_rho(
        &self,
        grids: Vec<[f64; 3]>,
        obt: u32,
        level: u32,
    ) -> Vec<(f64, [f64; 3], [[f64; 3]; 3])> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.obt_rho(&grid, obt, level))
            .collect()
    }

    pub fn RDG(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.RDG(&grid)).collect()
    }

    pub fn IRI(&self, grids: Vec<[f64; 3]>, alp: Option<f64>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.IRI(&grid, alp))
            .collect()
    }

    pub fn laplace(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.laplace(&grid))
            .collect()
    }

    pub fn signl2rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.signl2rho(&grid))
            .collect()
    }

    pub fn del_rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.del_rho(&grid))
            .collect()
    }

    pub fn fukui(&self, grids: Vec<[f64; 3]>, moleN: Mole, moleP: Mole) -> Vec<[f64; 4]> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.fukui(&grid, &moleN.inner, &moleP.inner))
            .collect()
    }

    pub fn lagkin(&self, grids: Vec<[f64; 3]>) -> Vec<[f64; 4]> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.lagkin(&grid))
            .collect()
    }

    pub fn hamkin(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.hamkin(&grid))
            .collect()
    }

    pub fn KED(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.KED(&grid)).collect()
    }

    pub fn ELF(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.ELF(&grid)).collect()
    }

    pub fn LOL(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.LOL(&grid)).collect()
    }

    pub fn infoentro(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.infoentro(&grid))
            .collect()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "density")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
