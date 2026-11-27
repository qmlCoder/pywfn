use crate::base::Mole;
use crate::props::gridprop::vals2numpy;

use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::{prelude::*, BoundObject};
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
    /// 创建一个新的密度计算器。
    ///
    /// * `mole` - 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn mol_rho_cm(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.mol_rho_cm(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    /// 计算开壳层（自旋非限制）分子在给定格点上的电子密度。
    ///
    /// * `py` - Python 解释器标记
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    /// * `level` - 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）
    ///
    /// 返回 `(密度值, 梯度, Hessian)` 的 NumPy 数组。
    pub fn mol_rho_dm(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.mol_rho_dm(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    pub fn pi_rho_pocv(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.pi_rho_pocv(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    pub fn pi_rho_mocv(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.pi_rho_mocv(&grid, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    pub fn ato_rho(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        level: usize,
    ) -> Vec<(Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>)> {
        let caler = self.caler();
        let ngrid = grids.len();
        let nato = self.mole.basis().atos().len("raw");
        let mut vals = vec![vec![rswfn::props::gridprop::GridVals::non(); ngrid]; nato];
        for (g, grid) in grids.iter().enumerate() {
            let ato_vals = caler.ato_rho(grid, level);
            for (a, ato_val) in ato_vals.into_iter().enumerate() {
                vals[a][g] = ato_val.clone();
            }
        }
        vals.iter().map(|val| vals2numpy(py, val, level)).collect()
    }

    /// 计算投影电子密度（仅含部分原子轨道贡献）。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的投影密度值。
    pub fn pro_rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.pro_rho(&grid))
            .collect()
    }

    /// 计算指定分子轨道在给定格点上的电子密度。
    ///
    /// * `py` - Python 解释器标记
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    /// * `obt` - 分子轨道索引（从 0 开始）
    /// * `level` - 计算级别（0=仅密度，1=密度+梯度，2=密度+梯度+Hessian）
    ///
    /// 返回 `(密度值, 梯度, Hessian)` 的 NumPy 数组。
    pub fn obt_rho(
        &self,
        py: Python,
        grids: Vec<[f64; 3]>,
        obt: usize,
        level: usize,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
        let caler = self.caler();
        let vals = grids
            .into_par_iter()
            .map(|grid| caler.obt_rho(&grid, obt, level))
            .collect();
        vals2numpy(py, &vals, level)
    }

    /// 计算约化密度梯度（Reduced Density Gradient）。
    ///
    /// 用于描述弱相互作用区域，是 RDG / NCI 分析的核心量。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 RDG 值。
    pub fn RDG(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.RDG(&grid)).collect()
    }

    /// 计算 IRI（Interaction Region Indicator）值。
    ///
    /// IRI 用于可视化分子内和分子间的弱相互作用区域。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    /// * `alp` - 可选参数，指数缩放因子；若为 `None` 则使用默认值
    ///
    /// 返回每个格点上的 IRI 值。
    pub fn IRI(&self, grids: Vec<[f64; 3]>, alp: Option<f64>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.IRI(&grid, alp))
            .collect()
    }

    /// 计算电子密度的拉普拉斯算子。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的拉普拉斯值。
    pub fn laplace(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.laplace(&grid))
            .collect()
    }

    /// 计算 sign(lambda2)*rho 值。
    ///
    /// 该量用于 NCI（非共价相互作用）分析，符号由 Hessian 矩阵的
    /// 第二特征值 lambda2 决定。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 sign(lambda2)*rho 值。
    pub fn signl2rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.signl2rho(&grid))
            .collect()
    }

    /// 计算电子密度的梯度模长 ||∇ρ||。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的梯度模长。
    pub fn del_rho(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.del_rho(&grid))
            .collect()
    }

    /// 计算 Fukui 函数值（亲核、亲电、自由基攻击）。
    ///
    /// 利用 N-1 和 N+1 电子态的分子计算 Fukui 函数。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    /// * `moleN` - N-1 电子态的分子对象
    /// * `moleP` - N+1 电子态的分子对象
    ///
    /// 返回每个格点上的 `[f+, f-, f0, ...]` 值数组。
    pub fn fukui(&self, grids: Vec<[f64; 3]>, moleN: Mole, moleP: Mole) -> Vec<[f64; 4]> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.fukui(&grid, &moleN.core, &moleP.core))
            .collect()
    }

    /// 计算 Lagrangian 动能密度。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 `[K, Kxx, Kyy, Kzz]` 值。
    pub fn lagkin(&self, grids: Vec<[f64; 3]>) -> Vec<[f64; 4]> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.lagkin(&grid))
            .collect()
    }

    /// 计算 Hamiltonian 动能密度。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 Hamiltonian 动能密度值。
    pub fn hamkin(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.hamkin(&grid))
            .collect()
    }

    /// 计算动能密度（Kinetic Energy Density）。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的动能密度值。
    pub fn KED(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.KED(&grid)).collect()
    }

    /// 计算电子局域化函数（Electron Localization Function）。
    ///
    /// ELF 用于描述电子的局域化程度，取值范围通常为 [0, 1]。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 ELF 值。
    pub fn ELF(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.ELF(&grid)).collect()
    }

    /// 计算局域化轨道指示函数（Localized Orbital Locator）。
    ///
    /// LOL 与 ELF 类似，用于表征电子的局域化特征。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的 LOL 值。
    pub fn LOL(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids.into_par_iter().map(|grid| caler.LOL(&grid)).collect()
    }

    /// 计算信息熵密度（Information Entropy Density）。
    ///
    /// * `grids` - 三维坐标点列表 `[[x, y, z], ...]`
    ///
    /// 返回每个格点上的信息熵密度值。
    pub fn infoentro(&self, grids: Vec<[f64; 3]>) -> Vec<f64> {
        let caler = self.caler();
        grids
            .into_par_iter()
            .map(|grid| caler.infoentro(&grid))
            .collect()
    }
    /// 获取 π 电子（POCV 方法）对应的分子对象
    pub fn mole_pi_pocv(&self) -> Mole {
        Mole {
            core: self.caler().mole_pi_pocv().clone(),
        }
    }

    /// 获取 π 电子（MOCV 方法）对应的分子对象
    pub fn mole_pi_mocv(&self) -> Mole {
        Mole {
            core: self.caler().mole_pi_mocv().clone(),
        }
    }

    /// 获取自旋密度对应的两个分子对象（alpha, beta）
    pub fn mole_spin(&self) -> (Mole, Mole) {
        let caler = self.caler();
        let (a, b) = caler.mole_spin();
        (Mole { core: a.clone() }, Mole { core: b.clone() })
    }

    /// 计算自旋密度
    ///
    /// ## 参数
    /// - `grid`: 格点坐标
    /// - `level`: 计算级别
    ///
    /// ## 返回
    /// 返回 (alpha 密度, beta 密度)
    pub fn mol_rho_spin(
        &self,
        grid: [f64; 3],
        level: usize,
    ) -> (
        crate::props::gridprop::GridVals,
        crate::props::gridprop::GridVals,
    ) {
        let [a, b] = self.caler().mol_rho_spin(&grid, level);
        (
            crate::props::gridprop::wrap_gridvals(a),
            crate::props::gridprop::wrap_gridvals(b),
        )
    }

    /// 计算指定格点上的活化能指标
    pub fn activate(&self, grid: [f64; 3]) -> f64 {
        self.caler().activate(&grid)
    }

    /// 计算指定格点上的曲率
    pub fn curvature(&self, grid: [f64; 3]) -> f64 {
        self.caler().curvature(&grid)
    }

    /// 计算指定格点上的核静电势
    pub fn nuc_esp(&self, grid: [f64; 3]) -> f64 {
        self.caler().nuc_esp(&grid)
    }

    /// 计算指定格点上的电子静电势
    pub fn ele_esp(&self, grid: [f64; 3]) -> f64 {
        self.caler().ele_esp(&grid)
    }

    /// 计算指定格点上的总静电势
    pub fn tot_esp(&self, grid: [f64; 3]) -> f64 {
        self.caler().tot_esp(&grid)
    }
}

/// 由原子波函数数据计算分子电子密度
///
/// ## 参数
/// - to_atms: 每个基函数所属的原子索引
/// - dmat_car: 笛卡尔基下的密度矩阵
/// - to_wfn0: 每个基函数在格点上的值
/// - to_wfn1: 每个基函数在格点上的梯度
/// - to_wfn2: 每个基函数在格点上的 Hessian
/// - tms: 需要参与计算的原子索引
/// - level: 计算级别
///
/// ## 返回
/// 返回 (密度值, 梯度, Hessian)
#[pyfunction]
pub fn mol_rho_dm(
    ato_atms: Vec<usize>,
    dmat_car: Bound<'_, PyArray2<f64>>,
    ato_wfn0: Vec<f64>,
    ato_wfn1: Vec<[f64; 3]>,
    ato_wfn2: Vec<[[f64; 3]; 3]>,
    atms: Vec<usize>,
    level: usize,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let dmat_car = dmat_car.readonly().as_array().to_owned();
    rswfn::props::gridprop::density::mol_rho_dm(
        &ato_atms, &dmat_car, &ato_wfn0, &ato_wfn1, &ato_wfn2, &atms, level,
    )
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "density")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(mol_rho_dm, &m)?)?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod density {
    #[pymodule_export]
    use super::Calculator;
}
