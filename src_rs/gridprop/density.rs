use crate::base::Mole;
use ndarray::{Array2, Array3};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3};
use pyo3::{prelude::*, BoundObject};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::gridprop::density::Calculator<'_> {
        rswfn::gridprop::density::Calculator::new(&self.mole.core)
    }
}

fn vals_split(
    vals: &Vec<rswfn::gridprop::GridVals>,
    level: usize,
) -> (Vec<f64>, Vec<Grad>, Vec<Hess>) {
    let mut vals0 = vec![];
    let mut vals1 = vec![];
    let mut vals2 = vec![];
    for val in vals {
        vals0.push(val.val0);
        if level > 0 {
            vals1.push(val.val1);
        }
        if level > 1 {
            vals2.push(val.val2);
        }
    }
    (vals0, vals1, vals2)
}

type Grad = [f64; 3];
type Hess = [[f64; 3]; 3];
fn grads2array(grads: Vec<Grad>) -> Array2<f64> {
    let ngrad = grads.len();
    let datas = grads.into_flattened();
    Array2::from_shape_vec((ngrad, 3), datas).unwrap()
}
fn hesss2array(hesss: Vec<Hess>) -> Array3<f64> {
    let nhess = hesss.len();
    let datas = hesss.into_flattened().into_flattened();
    Array3::from_shape_vec((nhess, 3, 3), datas).unwrap()
}

fn vals2numpy(
    py: Python,
    vals: &Vec<rswfn::gridprop::GridVals>,
    level: usize,
) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
    let (vals, grad, hess) = vals_split(vals, level);
    let vals = vals.into_pyarray(py).unbind();
    let grad = grads2array(grad).into_pyarray(py).unbind();
    let hess = hesss2array(hess).into_pyarray(py).unbind();
    (vals, grad, hess)
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mol_rho_cm(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.mol_rho_cm(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    pub fn mol_rho_dm(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.mol_rho_dm(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    // 获取原子轨道的格点
    pub fn ato_rho(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> Vec<(Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>)> {
        let caler = self.caler();
        let ngrid = grids.len();
        let nato = self.mole.basis().atos().len("raw");
        let mut vals = vec![vec![rswfn::gridprop::GridVals::non(); ngrid]; nato];
        for (g, grid) in grids.iter().enumerate() {
            let ato_vals = caler.ato_rho(grid, level);
            for (a, ato_val) in ato_vals.into_iter().enumerate() {
                vals[a][g] = ato_val.clone();
            }
        }
        vals.iter().map(|val| vals2numpy(py, val, level)).collect()
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
        py: Python,
        grids: Vec<[f64; 3]>,
        obt: usize,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.obt_rho(&grid, obt, level))
            .collect();
        vals2numpy(py, &vals, level)
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
            .map(|grid| caler.fukui(&grid, &moleN.core, &moleP.core))
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

#[pymodule]
pub mod density {
    #[pymodule_export]
    use super::Calculator;
}
