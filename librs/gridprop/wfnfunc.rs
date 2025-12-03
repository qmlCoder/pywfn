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
    fn caler(&self) -> rswfn::gridprop::wfnfunc::Calculator<'_> {
        rswfn::gridprop::wfnfunc::Calculator::new(&self.mole.inner)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn ato_wfn(
        &self,
        grids: Vec<[f64; 3]>, // 空间中任意一点的坐标
        level: u32,
        atms: Vec<u32>,
    ) -> Vec<(Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>)> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.ato_wfn(&grid, level, &atms))
            .collect()
    }

    // 计算一个点处的所有一个分子轨道波函数
    pub fn obt_wfn(
        &self,
        grids: Vec<[f64; 3]>, // 该点的坐标
        obt: u32,             // 分子轨道索引
        level: u32,
        atms: Vec<u32>,
    ) -> Vec<(f64, [f64; 3], [[f64; 3]; 3])> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| {
                let wfns = caler.ato_wfn(&grid, level, &atms);
                caler.obt_wfn(obt, level, &wfns)
            })
            .collect()
    }
}

#[pyfunction]
pub fn gtf(
    xyzs: Vec<[f64; 3]>, // 以原子为中心的格点坐标
    lmn: [u32; 3],
    alp: f64,
    level: u32,
) -> (Vec<f64>, Vec<[f64; 3]>, Vec<[[f64; 3]; 3]>) {
    let vals: Vec<(f64, [f64; 3], [[f64; 3]; 3])> = xyzs
        .into_par_iter()
        .map(|xyz| rswfn::gridprop::wfnfunc::gtf(&xyz, &lmn, alp, level))
        .collect();
    let val0: Vec<f64> = vals.iter().map(|val| val.0).collect();
    let mut val1: Vec<[f64; 3]> = Vec::new();
    let mut val2: Vec<[[f64; 3]; 3]> = Vec::new();
    if level >= 1 {
        val1 = vals.iter().map(|val| val.1).collect();
    }
    if level >= 2 {
        val2 = vals.iter().map(|val| val.2).collect();
    }
    (val0, val1, val2)
}

#[pyfunction]
// 计算收缩基函数（原子轨道），基函数的线性组合
pub fn cgf(
    xyzs: Vec<[f64; 3]>,
    lmn: [u32; 3],
    coes: Vec<f64>,
    alps: Vec<f64>,
    level: u32,
) -> Vec<(f64, [f64; 3], [[f64; 3]; 3])> {
    xyzs.into_par_iter()
        .map(|xyz| rswfn::gridprop::wfnfunc::cgf(&xyz, &lmn, &coes, &alps, level))
        .collect()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "wfnfunc")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(gtf, &m)?)?;
    parent_module.add_submodule(&m)
}
