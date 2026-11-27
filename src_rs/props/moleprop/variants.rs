use crate::base::Mole;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// 分子变体构造器。
///
/// 支持自旋分离（α/β）与 π 电子投影，生成新的分子对象。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::moleprop::variants::Calculator<'_> {
        let mut caler = rswfn::props::moleprop::variants::Calculator::new(&self.mole.core);
        caler.atms = &self.atms;
        caler
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的分子变体计算器。
    ///
    /// # 参数
    /// * `mole` - 分子对象
    /// * `atms` - 可选，参与 π 投影的原子索引列表，默认为全部原子
    #[new]
    pub fn new(mole: Mole, atms: Option<Vec<usize>>) -> Self {
        let atms = atms.unwrap_or_else(|| mole.core.atoms().idxs().clone());
        Self { mole, atms }
    }

    /// 将分子分离为 α 和 β 自旋两个独立的分子对象。
    pub fn mole_spin(&self) -> (Mole, Mole) {
        let (mole_a, mole_b) = self.caler().mole_spin();
        (Mole { core: mole_a }, Mole { core: mole_b })
    }

    /// 生成仅包含 π 电子成分（POCV 方法）的分子对象。
    pub fn mole_pi_pocv(&self) -> PyResult<Mole> {
        self.caler()
            .mole_pi_pocv()
            .map(|core| Mole { core })
            .map_err(|e| PyRuntimeError::new_err(e))
    }

    /// 生成仅包含 π 电子成分（MOCV 方法）的分子对象。
    pub fn mole_pi_mocv(&self) -> PyResult<Mole> {
        self.caler()
            .mole_pi_mocv()
            .map(|core| Mole { core })
            .map_err(|e| PyRuntimeError::new_err(e))
    }
}

/// 在父模块中注册 `variants` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "variants")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod variants {
    #[pymodule_export]
    use super::Calculator;
}
