use crate::base::Mole;
use ndarray::Array4;
use numpy::PyReadonlyArray4;
use pyo3::prelude::*;
use rswfn;

/// RHF（闭壳层 Hartree-Fock）计算器
///
/// 提供闭壳层自洽场（SCF）计算功能，用于求解分子的 Hartree-Fock 方程。
#[pyclass]
pub struct RHF {
    mole: Mole,
}

impl RHF {
    fn caler(&self) -> rswfn::qmcal::hf::RHF<'_> {
        rswfn::qmcal::hf::RHF::new(&self.mole.core)
    }
}

#[pymethods]
impl RHF {
    /// 创建一个新的 RHF 计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象。
    #[new]
    pub fn new(mole: Mole) -> PyResult<Self> {
        let rhf = Self { mole };
        Ok(rhf)
    }

    /// 执行 RHF 自洽场计算。
    ///
    /// 迭代求解 Hartree-Fock 方程，获得收敛的分子轨道和能量。
    /// 可选地提供双电子积分矩阵以减少重复计算。
    ///
    /// # 参数
    /// * `emat` - 可选，双电子积分矩阵（四维数组）。如果提供，则使用该矩阵
    ///            进行计算，避免重新计算双电子积分。
    pub fn run(&self, emat: Option<PyReadonlyArray4<f64>>) {
        let emat = match emat {
            Some(emat) => {
                let emat = emat.as_array().to_owned();
                let emat = rswfn::rsint::Mat2e::from_mat4(emat);
                Some(emat)
            }
            None => None,
        };
        self.caler().run(emat);
    }
}

/// 在父模块中注册 `hf` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "hf")?;
    m.add_class::<RHF>()?;
    parent_module.add_submodule(&m)
}
