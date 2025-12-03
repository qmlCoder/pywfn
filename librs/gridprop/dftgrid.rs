use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::gridprop::dftgrid::Calculator<'_> {
        rswfn::gridprop::dftgrid::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    // 计算dft原子径向格点
    pub fn rad_grids(
        &self,
        R: f64,      // 原子半径
        nrad: usize, // 径向格点数量
    ) -> (Vec<f64>, Vec<f64>) {
        self.caler().rad_grids(R, nrad)
    }

    pub fn atm_grids(
        &self,
        iatm: u32,
        nrad: usize, // 径向格点数量
        fsph: usize, // 球面格点函数
    ) -> (Vec<[f64; 3]>, Vec<f64>) {
        self.caler().atm_grids(iatm, nrad, fsph)
    }

    pub fn a2m_weits(
        &self,
        iatm: u32,                // 第多少个原子
        atm_grids: Vec<[f64; 3]>, // 原子格点
        atm_weits: Vec<f64>,      // 原子格点权重
    ) -> Vec<f64> {
        self.caler().a2m_weits(iatm, &atm_grids, &atm_weits)
    }

    pub fn mol_grids(&self, nrad: usize, fsph: usize) -> PyResult<(Vec<[f64; 3]>, Vec<f64>)> {
        let res = self.caler().mol_grids(nrad, fsph);
        Ok(res)
    }

    pub fn frag_grids(&self, frag: Vec<u32>) -> (Vec<[f64; 3]>, Vec<f64>) {
        self.caler().frag_grids(&frag)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "dftgrid")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
