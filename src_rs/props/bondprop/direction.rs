use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::bondprop::direction::Calculator<'_> {
        rswfn::props::bondprop::direction::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的键方向计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算两个原子间的垂直方向
    ///
    /// 返回两个原子之间垂直于键轴的方向向量，用于描述键的 π 方向。
    ///
    /// ## 参数
    /// - `ai`: 第一个原子的索引
    /// - `aj`: 第二个原子的索引
    ///
    /// ## 返回
    /// 垂直方向向量 [x, y, z] 的列表，通常包含两个正交方向
    pub fn verts(&self, ai: usize, aj: usize) -> Vec<[f64; 3]> {
        let dirs = self.caler().verts(ai, aj);
        dirs.into_iter().map(|dir| dir.into_arr()).collect()
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "direction")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

// #[pymodule]
// pub mod direction {
//     #[pymodule_export]
//     use super::Calculator;
// }
