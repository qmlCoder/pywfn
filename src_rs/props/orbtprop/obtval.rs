use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

/// 轨道数值计算器
///
/// 提供轨道相关数值量的计算功能，例如 π 电子投影值等。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::orbtprop::obtval::Calculator<'_> {
        rswfn::props::orbtprop::obtval::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的轨道数值计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算每个原子 π 电子投影系数（pocv 方法）。
    ///
    /// 使用轨道投影方法计算每个原子上 π 电子的分布情况。
    ///
    /// # 返回
    /// * `Vec<f64>` - 每个原子对应的 π 电子投影值列表。
    pub fn pi_eles_pocv(&self) -> Vec<f64> {
        self.caler().pi_eles_pocv()
    }

    /// 计算每个原子 π 电子投影系数（mocv 方法）。
    ///
    /// 使用分子轨道变换方法计算每个原子上 π 电子的分布情况。
    ///
    /// # 返回
    /// * `Vec<f64>` - 每个原子对应的 π 电子投影值列表。
    pub fn pi_eles_mocv(&self) -> Vec<f64> {
        self.caler().pi_eles_mocv()
    }
    /// 计算各原子轨道上的电子数
    ///
    /// ## 返回
    /// 每个原子的电子数列表
    pub fn eles(&self) -> Vec<f64> {
        self.caler().eles()
    }
}

/// 在父模块中注册 `obtval` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "obtval")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod obtval {
    #[pymodule_export]
    use super::Calculator;
}
