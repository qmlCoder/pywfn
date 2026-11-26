use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::fragprop::energy::Calculator<'_> {
        rswfn::props::fragprop::energy::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的碎片能量计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }
    /// 计算 EIEBA（Energy of Interaction between EDA and BAE）
    ///
    /// 计算两个碎片之间的相互作用能，包括静电、交换和轨道相互作用分量。
    ///
    /// ## 参数
    /// - `fragA`: 碎片 A 的原子索引列表
    /// - `fragB`: 碎片 B 的原子索引列表
    ///
    /// ## 返回
    /// - `f64`: 静电相互作用能
    /// - `f64`: 交换排斥能
    /// - `f64`: 轨道相互作用能
    pub fn EIEBA(&self, fragA: Vec<usize>, fragB: Vec<usize>) -> (f64, f64, f64) {
        self.caler().EIEBA(&fragA, &fragB)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "energy")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod energy {
    #[pymodule_export]
    use super::Calculator;
}
