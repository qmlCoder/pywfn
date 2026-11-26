use crate::base::Mole;
use pyo3::prelude::*;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use rswfn;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    /// 创建底层的 `rswfn::props::gridprop::density::Calculator` 实例。
    fn caler(&self) -> rswfn::props::gridprop::density::Calculator<'_> {
        rswfn::props::gridprop::density::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的势能计算器。
    ///
    /// * `mole` - 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算原子核势能（核-电子吸引势）在给定格点上的值。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的核势能值。
    pub fn nuc_potential(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.nuc_esp(&grid))
            .collect()
    }

    /// 计算电子势能（电子-电子排斥势）在给定格点上的值。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的电子势能值。
    pub fn ele_potential(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.ele_esp(&grid))
            .collect()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "potential")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod potential {
    #[pymodule_export]
    use super::Calculator;
}
