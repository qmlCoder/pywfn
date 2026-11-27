use crate::base::Mole;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// 键能量的可加分解结果（闭壳层，杂化泛函）。
///
/// 把完整电子能量按密度拆成"键"与"其余"两部分，满足
/// `E_bond + E_rest + E_cross = E_elec`。
#[pyclass]
#[derive(Clone, Copy)]
pub struct BondEnergy {
    core: rswfn::props::bondprop::energy::BondEnergy,
}

#[pymethods]
impl BondEnergy {
    /// 键自能量 `E_bond`
    #[getter]
    pub fn e_bond(&self) -> f64 {
        self.core.e_bond
    }
    /// 其余部分能量 `E_rest`
    #[getter]
    pub fn e_rest(&self) -> f64 {
        self.core.e_rest
    }
    /// 键-其余交叉项 `E_cross`
    #[getter]
    pub fn e_cross(&self) -> f64 {
        self.core.e_cross
    }
    /// 交叉项的库仑部分 `Tr(D_bond J_rest)`
    #[getter]
    pub fn e_cross_coul(&self) -> f64 {
        self.core.e_cross_coul
    }
    /// 交叉项的精确交换部分 `−aₓ/2·Tr(D_bond K_rest)`
    #[getter]
    pub fn e_cross_exch(&self) -> f64 {
        self.core.e_cross_exch
    }
    /// 交叉项的 XC 残差
    #[getter]
    pub fn e_cross_xc(&self) -> f64 {
        self.core.e_cross_xc
    }
    /// 完整电子能量 `E_elec`
    #[getter]
    pub fn e_elec(&self) -> f64 {
        self.core.e_elec
    }
    /// 核-核排斥能 `E_nuc`
    #[getter]
    pub fn e_nuc(&self) -> f64 {
        self.core.e_nuc
    }
    /// 总能量 `E_total = E_elec + E_nuc`
    #[getter]
    pub fn e_total(&self) -> f64 {
        self.core.e_total
    }

    /// 键对电子能量的可加贡献 `E_bond + E_cross`
    pub fn bond(&self) -> f64 {
        self.core.bond()
    }
}

/// 键能量计算器。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::bondprop::energy::Calculator<'_> {
        rswfn::props::bondprop::energy::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的键能量计算器。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 键（两原子）对电子能量的可加贡献 `E_bond + E_cross`。
    pub fn bond_energy(&self, atm1: usize, atm2: usize) -> PyResult<f64> {
        self.caler()
            .bond_energy(atm1, atm2)
            .map_err(|e| PyRuntimeError::new_err(e))
    }

    /// 对指定原子对的完整电子能量做 bond / rest / cross 分解。
    pub fn bond_parts(&self, atm1: usize, atm2: usize) -> PyResult<BondEnergy> {
        self.caler()
            .bond_parts(atm1, atm2)
            .map(|core| BondEnergy { core })
            .map_err(|e| PyRuntimeError::new_err(e))
    }
}

/// 在父模块中注册 `energy` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "energy")?;
    m.add_class::<Calculator>()?;
    m.add_class::<BondEnergy>()?;
    parent_module.add_submodule(&m)
}
