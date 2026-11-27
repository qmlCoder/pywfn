use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

/// 拓扑分析计算器
///
/// 从若干初猜点出发寻找分子电子密度（或 π 密度）的临界点，
/// 用于 AIM 等拓扑分析。
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::gridprop::topolog::Calculator<'_> {
        rswfn::props::gridprop::topolog::Calculator::new(&self.mole.core)
    }
}

/// 由密度函数名生成的格点函数
#[pyclass]
pub struct TopoFunc {
    pub mole: Mole,
    pub name: String,
}

#[pymethods]
impl TopoFunc {
    /// 计算指定格点上的函数值
    ///
    /// # 参数
    /// * `grid` - 格点坐标
    ///
    /// # 返回
    /// 格点上的函数值及其导数
    pub fn __call__(&self, grid: [f64; 3]) -> crate::props::gridprop::GridVals {
        let caler = rswfn::props::gridprop::topolog::Calculator::new(&self.mole.core);
        let func = caler.get_func(&self.name);
        crate::props::gridprop::wrap_gridvals(func(&grid))
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的拓扑分析计算器。
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 获取分子中的临界点初猜（原子中心 + 键中心）。
    pub fn get_grids(&self) -> Vec<[f64; 3]> {
        self.caler().get_grids()
    }

    /// 从指定格点出发寻找临界点
    ///
    /// # 参数
    /// * `grid` - 初始格点坐标
    /// * `func` - 密度函数名，可选 `"mol_rho_dm"` / `"pi_rho_pocv"` / `"pi_rho_mocv"`
    ///
    /// # 返回
    /// 返回 `(位置, 密度值, 临界点类型)`
    pub fn get_point(&self, grid: [f64; 3], func: &str) -> ([f64; 3], f64, i32) {
        let caler = self.caler();
        let f = caler.get_func(func);
        caler.get_point(&grid, &f)
    }

    /// 获取指定名称的密度函数对象
    ///
    /// # 参数
    /// * `func` - 密度函数名，可选 `"mol_rho_dm"` / `"pi_rho_pocv"` / `"pi_rho_mocv"`
    ///
    /// # 返回
    /// 可直接调用的函数对象，输入格点返回函数值
    pub fn get_func(&self, func: &str) -> TopoFunc {
        TopoFunc {
            mole: self.mole.clone(),
            name: func.to_string(),
        }
    }

    /// 执行拓扑分析，返回全部临界点。
    ///
    /// # 参数
    /// * `func` - 密度函数名，可选 `"mol_rho_dm"` / `"pi_rho_pocv"` / `"pi_rho_mocv"`
    ///
    /// # 返回
    /// 每个临界点的 `(位置 [x,y,z], 密度值, 临界点类型)` 列表；
    /// 类型为 Hessian 正、负特征值个数之差（3=核临界点，-3=笼临界点等）。
    pub fn analysis(&self, func: &str) -> Vec<([f64; 3], f64, i32)> {
        self.caler().analysis(func)
    }
}

/// 在父模块中注册 `topolog` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "topolog")?;
    m.add_class::<Calculator>()?;
    m.add_class::<TopoFunc>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod topolog {
    #[pymodule_export]
    use super::Calculator;
}
