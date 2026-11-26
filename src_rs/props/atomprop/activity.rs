use crate::base::{Mole, Stm};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn::{self, base::Vec3};
use std::collections::HashMap;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::atomprop::activity::Calculator<'_> {
        rswfn::props::atomprop::activity::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的活性计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn fukui(&self, mole_n: Mole, mole_p: Mole, ctype: &str) -> PyResult<Vec<[f64; 7]>> {
        let vals = self.caler().fukui(&mole_n.core, &mole_p.core, ctype);
        Ok(vals)
    }

    pub fn fukui_pi_pocv(
        &self,
        mole_n: Mole,
        mole_p: Mole,
        ctype: &str,
    ) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<[f64; 7]>)> {
        let (dirs, vals) = self
            .caler()
            .fukui_pi_pocv(&mole_n.core, &mole_p.core, ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let dirs = dirs
            .iter()
            .map(|(key, dir)| (*key, dir.into_arr()))
            .collect();
        Ok((dirs, vals))
    }

    pub fn fukui_pi_mocv(
        &self,
        mole_n: Mole,
        mole_p: Mole,
        ctype: &str,
    ) -> PyResult<(HashMap<usize, Stm>, Vec<usize>, Vec<[f64; 7]>)> {
        let (stms, atos, vals) = self
            .caler()
            .fukui_pi_mocv(&mole_n.core, &mole_p.core, ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let stms = stms
            .iter()
            .map(|(key, stm)| (*key, Stm { core: stm.clone() }))
            .collect();
        Ok((stms, atos, vals))
    }

    pub fn fukui_dir(
        &self,
        atm: usize,
        dir: [f64; 3],
        mole_n: Mole,
        mole_p: Mole,
        ctype: &str,
    ) -> PyResult<[f64; 7]> {
        let dir = Vec3::from_arr(dir);
        let vals = self
            .caler()
            .fukui_dir(atm, &dir, &mole_n.core, &mole_p.core, ctype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        Ok(vals)
    }

    pub fn freev(&self, otype: &str) -> Vec<[f64; 3]> {
        self.caler().freev(otype)
    }

    pub fn freev_pi(&self, otype: &str) -> PyResult<(HashMap<usize, [f64; 3]>, Vec<[f64; 3]>)> {
        let (dirs, vals) = self
            .caler()
            .freev_pi(otype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))?;
        let dirs = dirs
            .iter()
            .map(|(atm, dir)| (*atm, dir.into_arr()))
            .collect();
        Ok((dirs, vals))
    }

    pub fn freev_dir(
        &self,
        atm: usize,
        dir: [f64; 3],
        nebs: Vec<usize>,
        otype: &str,
        mtype: &str,
    ) -> PyResult<(f64, Vec<f64>, f64)> {
        let dir = Vec3::from_arr(dir);
        // rswfn 已将该方法拆分为 pocv/mocv，默认沿用旧实现的 pocv
        self.caler()
            .freev_dir(atm, &dir, &nebs, otype, mtype)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "activity")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
