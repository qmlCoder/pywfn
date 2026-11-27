use crate::base::{Basis, Coefs, Mole};
use ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// 把格点坐标列表转换为 `(ngrid, 3)` 的 NumPy 数组
pub(crate) fn grids2array(py: Python, grids: &[[f64; 3]]) -> Py<PyArray2<f64>> {
    let ngrid = grids.len();
    let flat: Vec<f64> = grids.iter().flat_map(|grid| grid.iter().copied()).collect();
    Array2::from_shape_vec((ngrid, 3), flat)
        .expect("格点坐标数组形状不匹配")
        .into_pyarray(py)
        .unbind()
}

/// 闭壳层 Kohn-Sham（DFT）计算结果
#[pyclass]
pub struct RKSResult {
    core: rswfn::qmcal::ks::RKSResult,
}

#[pymethods]
impl RKSResult {
    /// 电子能量（含交换相关能）
    #[getter]
    pub fn ele_eng(&self) -> f64 {
        self.core.ele_eng
    }
    /// 交换相关能
    #[getter]
    pub fn xc_eng(&self) -> f64 {
        self.core.xc_eng
    }
    /// 核-核排斥能
    #[getter]
    pub fn nuc_eng(&self) -> f64 {
        self.core.nuc_eng
    }
    /// 泛函名称
    #[getter]
    pub fn func(&self) -> String {
        self.core.func.name().to_string()
    }
    /// 是否自洽收敛
    #[getter]
    pub fn converged(&self) -> bool {
        self.core.converged
    }
    /// 实际迭代次数
    #[getter]
    pub fn niter(&self) -> usize {
        self.core.niter
    }
    /// 最后一次迭代的能量变化（Hartree）
    #[getter]
    pub fn delta_eng(&self) -> f64 {
        self.core.delta_eng
    }
    /// 分子轨道系数
    #[getter]
    pub fn coefs(&self) -> Coefs {
        Coefs {
            core: self.core.coefs.clone(),
        }
    }
    /// 基组信息
    #[getter]
    pub fn basis(&self) -> Basis {
        Basis {
            core: self.core.basis.clone(),
        }
    }
    /// 实际使用的积分形式（"sph" 或 "car"）
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }
    /// 数值积分格点数量
    #[getter]
    pub fn ngrid(&self) -> usize {
        self.core.ngrid
    }
    /// 单电子积分耗时（秒）
    #[getter]
    pub fn time_int1e(&self) -> f64 {
        self.core.time_int1e.as_secs_f64()
    }
    /// 双电子积分耗时（秒）
    #[getter]
    pub fn time_int2e(&self) -> f64 {
        self.core.time_int2e.as_secs_f64()
    }
    /// 自洽场迭代总耗时（秒）
    #[getter]
    pub fn time_scf(&self) -> f64 {
        self.core.time_scf.as_secs_f64()
    }
    /// 构建 J/K 矩阵耗时（秒）
    #[getter]
    pub fn time_jk(&self) -> f64 {
        self.core.time_jk.as_secs_f64()
    }
    /// 数值积分耗时（秒）
    #[getter]
    pub fn time_grid(&self) -> f64 {
        self.core.time_grid.as_secs_f64()
    }
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.core.total_eng()
    }
}

/// 自旋极化（开壳层）Kohn-Sham（DFT）计算结果
#[pyclass]
pub struct UKSResult {
    core: rswfn::qmcal::ks::UKSResult,
}

#[pymethods]
impl UKSResult {
    #[getter]
    pub fn ele_eng(&self) -> f64 {
        self.core.ele_eng
    }
    #[getter]
    pub fn xc_eng(&self) -> f64 {
        self.core.xc_eng
    }
    #[getter]
    pub fn nuc_eng(&self) -> f64 {
        self.core.nuc_eng
    }
    /// 泛函名称
    #[getter]
    pub fn func(&self) -> String {
        self.core.func.name().to_string()
    }
    #[getter]
    pub fn converged(&self) -> bool {
        self.core.converged
    }
    #[getter]
    pub fn niter(&self) -> usize {
        self.core.niter
    }
    #[getter]
    pub fn delta_eng(&self) -> f64 {
        self.core.delta_eng
    }
    /// 分子轨道系数
    #[getter]
    pub fn coefs(&self) -> Coefs {
        Coefs {
            core: self.core.coefs.clone(),
        }
    }
    /// 基组信息
    #[getter]
    pub fn basis(&self) -> Basis {
        Basis {
            core: self.core.basis.clone(),
        }
    }
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }
    #[getter]
    pub fn ngrid(&self) -> usize {
        self.core.ngrid
    }
    /// α 电子数
    #[getter]
    pub fn nela(&self) -> usize {
        self.core.nela
    }
    /// β 电子数
    #[getter]
    pub fn nelb(&self) -> usize {
        self.core.nelb
    }
    /// ⟨S²⟩ 期望值
    #[getter]
    pub fn spin_sq(&self) -> f64 {
        self.core.spin_sq
    }
    /// 单电子积分耗时（秒）
    #[getter]
    pub fn time_int1e(&self) -> f64 {
        self.core.time_int1e.as_secs_f64()
    }
    /// 双电子积分耗时（秒）
    #[getter]
    pub fn time_int2e(&self) -> f64 {
        self.core.time_int2e.as_secs_f64()
    }
    /// 自洽场迭代总耗时（秒）
    #[getter]
    pub fn time_scf(&self) -> f64 {
        self.core.time_scf.as_secs_f64()
    }
    /// 构建 J/K 矩阵耗时（秒）
    #[getter]
    pub fn time_jk(&self) -> f64 {
        self.core.time_jk.as_secs_f64()
    }
    /// 数值积分耗时（秒）
    #[getter]
    pub fn time_grid(&self) -> f64 {
        self.core.time_grid.as_secs_f64()
    }
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.core.total_eng()
    }
    /// 自旋多重度 `2S + 1`
    pub fn multip(&self) -> usize {
        self.core.multip()
    }
}

/// 闭壳层 Kohn-Sham（DFT）计算器
#[pyclass]
pub struct RKS {
    pub mole: Mole,
    pub func: String,
    pub grid: Option<(usize, usize)>,
    /// 是否对格点做裁剪；`None` 表示沿用上游默认行为
    pub prune: Option<bool>,
    pub opt_int_form: Option<String>,
    pub opt_max_iter: Option<usize>,
    pub opt_tol: Option<f64>,
    pub opt_screen: Option<f64>,
    pub opt_max_mem: Option<usize>,
    pub opt_ao_cut: Option<f64>,
}

impl RKS {
    /// 用统一入口构造对象，避免各处结构体字面量遗漏字段
    fn build(mole: Mole, func: String, grid: Option<(usize, usize)>, prune: Option<bool>) -> Self {
        Self {
            mole,
            func,
            grid,
            prune,
            opt_int_form: None,
            opt_max_iter: None,
            opt_tol: None,
            opt_screen: None,
            opt_max_mem: None,
            opt_ao_cut: None,
        }
    }

    fn caler(&self) -> rswfn::qmcal::ks::RKS<'_> {
        let func = rswfn::rsxcs::XcFunc::from_name(&self.func).expect("未知泛函");
        let mut caler = match self.grid {
            None => rswfn::qmcal::ks::RKS::new(&self.mole.core, func),
            Some((radial, angular)) => match self.prune {
                None => rswfn::qmcal::ks::RKS::with_grid(&self.mole.core, radial, angular, func),
                Some(prune) => rswfn::qmcal::ks::RKS::with_grid_prune(
                    &self.mole.core,
                    radial,
                    angular,
                    func,
                    prune,
                ),
            },
        };
        if let Some(value) = &self.opt_int_form {
            caler.int_form = value.clone();
        }
        if let Some(value) = self.opt_max_iter {
            caler.max_iter = value;
        }
        if let Some(value) = self.opt_tol {
            caler.tol = value;
        }
        if let Some(value) = self.opt_screen {
            caler.screen = value;
        }
        if let Some(value) = self.opt_max_mem {
            caler.max_mem = value;
        }
        if let Some(value) = self.opt_ao_cut {
            caler.ao_cut = value;
        }
        caler
    }

    fn check(func: &str) -> PyResult<()> {
        if rswfn::rsxcs::XcFunc::from_name(func).is_none() {
            return Err(PyRuntimeError::new_err(format!("未知泛函: {func}")));
        }
        Ok(())
    }
}

#[pymethods]
impl RKS {
    /// 按泛函名创建闭壳层 KS 计算器（默认格点）。
    #[new]
    pub fn new(mole: Mole, func: String) -> PyResult<Self> {
        Self::check(&func)?;
        Ok(Self::build(mole, func, None, None))
    }

    /// 按泛函名创建闭壳层 KS 计算器（默认格点）；泛函名无效时返回 `None`。
    #[staticmethod]
    pub fn from_name(mole: Mole, name: String) -> Option<Self> {
        rswfn::rsxcs::XcFunc::from_name(&name)?;
        Some(Self::build(mole, name, None, None))
    }

    /// 用指定的径向/球面格点数创建计算器。
    #[staticmethod]
    pub fn with_grid(mole: Mole, radial: usize, angular: usize, func: String) -> PyResult<Self> {
        Self::check(&func)?;
        Ok(Self::build(mole, func, Some((radial, angular)), None))
    }

    /// 用指定的径向/球面格点数创建计算器，并指定是否裁剪格点。
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `radial`: 径向格点数
    /// - `angular`: 球面格点数
    /// - `func`: 泛函名（如 "b3lyp"）
    /// - `prune`: 是否对格点做裁剪
    #[staticmethod]
    pub fn with_grid_prune(
        mole: Mole,
        radial: usize,
        angular: usize,
        func: String,
        prune: bool,
    ) -> PyResult<Self> {
        Self::check(&func)?;
        Ok(Self::build(mole, func, Some((radial, angular)), Some(prune)))
    }

    /// 执行自洽场计算。
    pub fn run(&self) -> RKSResult {
        RKSResult {
            core: self.caler().run(),
        }
    }

    /// 核-核排斥能。
    pub fn nuc_eng(&self) -> f64 {
        self.caler().nuc_eng()
    }

    /// 原子轨道统计信息：`(屏幕后有效原子轨道数, 最大系数绝对值, 原始原子轨道数)`
    pub fn ao_stats(&self) -> (usize, f64, usize) {
        self.caler().ao_stats()
    }

    /// 在数值积分格点上计算电子密度及其梯度
    ///
    /// ## 参数
    /// - `dmat`: 密度矩阵（NumPy 二维数组）
    ///
    /// ## 返回
    /// 三元组 `(格点密度, 格点梯度, 格点密度平方模)`，梯度为 `(ngrid, 3)` 数组
    pub fn calc_rho(
        &self,
        py: Python,
        dmat: Bound<'_, PyArray2<f64>>,
    ) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray1<f64>>) {
        let dmat = dmat.readonly().as_array().to_owned();
        let caler = self.caler();
        let (rho, grad, sigma) = caler.calc_rho(&dmat);
        (
            rho.into_pyarray(py).unbind(),
            grids2array(py, &grad),
            sigma.into_pyarray(py).unbind(),
        )
    }

    /// 计算交换相关能
    ///
    /// ## 参数
    /// - `dmat`: 密度矩阵（NumPy 二维数组）
    ///
    /// ## 返回
    /// 交换相关能（Hartree）
    pub fn xc_eng(&self, dmat: Bound<'_, PyArray2<f64>>) -> f64 {
        let dmat = dmat.readonly().as_array().to_owned();
        self.caler().xc_eng(&dmat)
    }

    /// 计算电子在数值积分格点上的静电相互作用能
    ///
    /// ## 参数
    /// - `dmat`: 密度矩阵（NumPy 二维数组）
    ///
    /// ## 返回
    /// 格点贡献的电子能量（Hartree）
    pub fn grid_ele(&self, dmat: Bound<'_, PyArray2<f64>>) -> f64 {
        let dmat = dmat.readonly().as_array().to_owned();
        self.caler().grid_ele(&dmat)
    }

    /// 泛函名称
    #[getter]
    pub fn func(&self) -> String {
        self.func.clone()
    }

    /// 积分形式（"sph" 或 "car"）
    #[getter]
    pub fn int_form(&self) -> String {
        self.caler().int_form.clone()
    }
    #[setter]
    pub fn set_int_form(&mut self, value: String) {
        self.opt_int_form = Some(value);
    }

    /// 径向格点数
    #[getter]
    pub fn radial(&self) -> usize {
        self.caler().radial
    }

    /// 球面格点数
    #[getter]
    pub fn angular(&self) -> usize {
        self.caler().angular
    }

    /// 最大 SCF 迭代次数
    #[getter]
    pub fn max_iter(&self) -> usize {
        self.caler().max_iter
    }
    #[setter]
    pub fn set_max_iter(&mut self, value: usize) {
        self.opt_max_iter = Some(value);
    }

    /// SCF 收敛阈值
    #[getter]
    pub fn tol(&self) -> f64 {
        self.caler().tol
    }
    #[setter]
    pub fn set_tol(&mut self, value: f64) {
        self.opt_tol = Some(value);
    }

    /// 数值积分格点坐标 `(ngrid, 3)`
    #[getter]
    pub fn grids(&self, py: Python) -> Py<PyArray2<f64>> {
        let caler = self.caler();
        grids2array(py, &caler.grids)
    }

    /// 数值积分格点权重
    #[getter]
    pub fn weits(&self, py: Python) -> Py<PyArray1<f64>> {
        self.caler().weits.clone().into_pyarray(py).unbind()
    }

    /// 数值积分格点数量
    #[getter]
    pub fn ngrid(&self) -> usize {
        self.caler().ngrid
    }

    /// 原子轨道数量
    #[getter]
    pub fn nato(&self) -> usize {
        self.caler().nato
    }

    /// 原子轨道在格点上的取值
    #[getter]
    pub fn chi0(&self, py: Python) -> Py<PyArray1<f64>> {
        self.caler().chi0.clone().into_pyarray(py).unbind()
    }

    /// 原子轨道在格点上的梯度分量
    #[getter]
    pub fn chi1(&self, py: Python) -> Py<PyArray1<f64>> {
        self.caler().chi1.clone().into_pyarray(py).unbind()
    }

    /// 双电子积分 screening 阈值
    #[getter]
    pub fn screen(&self) -> f64 {
        self.caler().screen
    }
    #[setter]
    pub fn set_screen(&mut self, value: f64) {
        self.opt_screen = Some(value);
    }

    /// 内存上限（字节）
    #[getter]
    pub fn max_mem(&self) -> usize {
        self.caler().max_mem
    }
    #[setter]
    pub fn set_max_mem(&mut self, value: usize) {
        self.opt_max_mem = Some(value);
    }

    /// 原子轨道在格点上的截断阈值
    #[getter]
    pub fn ao_cut(&self) -> f64 {
        self.caler().ao_cut
    }
    #[setter]
    pub fn set_ao_cut(&mut self, value: f64) {
        self.opt_ao_cut = Some(value);
    }
}

/// 自旋极化（开壳层）Kohn-Sham（DFT）计算器
#[pyclass]
pub struct UKS {
    pub mole: Mole,
    pub func: String,
    pub grid: Option<(usize, usize)>,
    /// 是否对格点做裁剪；`None` 表示沿用上游默认行为
    pub prune: Option<bool>,
}

impl UKS {
    fn build(mole: Mole, func: String, grid: Option<(usize, usize)>, prune: Option<bool>) -> Self {
        Self {
            mole,
            func,
            grid,
            prune,
        }
    }

    fn caler(&self) -> rswfn::qmcal::ks::UKS<'_> {
        let func = rswfn::rsxcs::XcFunc::from_name(&self.func).expect("未知泛函");
        match self.grid {
            None => rswfn::qmcal::ks::UKS::new(&self.mole.core, func),
            Some((radial, angular)) => match self.prune {
                None => rswfn::qmcal::ks::UKS::with_grid(&self.mole.core, radial, angular, func),
                Some(prune) => rswfn::qmcal::ks::UKS::with_grid_prune(
                    &self.mole.core,
                    radial,
                    angular,
                    func,
                    prune,
                ),
            },
        }
    }
}

#[pymethods]
impl UKS {
    /// 按泛函名创建开壳层 KS 计算器（默认格点）。
    #[new]
    pub fn new(mole: Mole, func: String) -> PyResult<Self> {
        if rswfn::rsxcs::XcFunc::from_name(&func).is_none() {
            return Err(PyRuntimeError::new_err(format!("未知泛函: {func}")));
        }
        Ok(Self::build(mole, func, None, None))
    }

    /// 按泛函名创建开壳层 KS 计算器（默认格点）；泛函名无效时返回 `None`。
    #[staticmethod]
    pub fn from_name(mole: Mole, name: String) -> Option<Self> {
        rswfn::rsxcs::XcFunc::from_name(&name)?;
        Some(Self::build(mole, name, None, None))
    }

    /// 用指定的径向/球面格点数创建计算器。
    #[staticmethod]
    pub fn with_grid(mole: Mole, radial: usize, angular: usize, func: String) -> PyResult<Self> {
        if rswfn::rsxcs::XcFunc::from_name(&func).is_none() {
            return Err(PyRuntimeError::new_err(format!("未知泛函: {func}")));
        }
        Ok(Self::build(mole, func, Some((radial, angular)), None))
    }

    /// 用指定的径向/球面格点数创建计算器，并指定是否裁剪格点。
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `radial`: 径向格点数
    /// - `angular`: 球面格点数
    /// - `func`: 泛函名（如 "b3lyp"）
    /// - `prune`: 是否对格点做裁剪
    #[staticmethod]
    pub fn with_grid_prune(
        mole: Mole,
        radial: usize,
        angular: usize,
        func: String,
        prune: bool,
    ) -> PyResult<Self> {
        if rswfn::rsxcs::XcFunc::from_name(&func).is_none() {
            return Err(PyRuntimeError::new_err(format!("未知泛函: {func}")));
        }
        Ok(Self::build(mole, func, Some((radial, angular)), Some(prune)))
    }

    /// 执行自洽场计算。
    pub fn run(&self) -> UKSResult {
        UKSResult {
            core: self.caler().run(),
        }
    }

    /// 核-核排斥能。
    pub fn nuc_eng(&self) -> f64 {
        self.caler().nuc_eng()
    }

    /// 在数值积分格点上计算 α/β 自旋密度及其梯度
    ///
    /// ## 参数
    /// - `dmat_a`: α 密度矩阵（NumPy 二维数组）
    /// - `dmat_b`: β 密度矩阵（NumPy 二维数组）
    ///
    /// ## 返回
    /// 自旋密度对象（含密度、梯度及各密度梯度内积）
    pub fn calc_rho_spin(
        &self,
        dmat_a: Bound<'_, PyArray2<f64>>,
        dmat_b: Bound<'_, PyArray2<f64>>,
    ) -> SpinRho {
        let dmat_a = dmat_a.readonly().as_array().to_owned();
        let dmat_b = dmat_b.readonly().as_array().to_owned();
        SpinRho {
            core: self.caler().calc_rho_spin(&dmat_a, &dmat_b),
        }
    }

    /// 底层闭壳层 KS 计算器（对应上游字段 `base`）
    #[getter]
    pub fn base(&self) -> RKS {
        RKS::build(self.mole.clone(), self.func.clone(), self.grid, self.prune)
    }

    /// α 电子数
    #[getter]
    pub fn nela(&self) -> usize {
        self.caler().nela
    }

    /// β 电子数
    #[getter]
    pub fn nelb(&self) -> usize {
        self.caler().nelb
    }
}

/// 数值积分格点分块信息
#[pyclass]
pub struct GridBlock {
    core: rswfn::qmcal::ks::GridBlock,
}

#[pymethods]
impl GridBlock {
    /// 创建一个格点分块
    ///
    /// ## 参数
    /// - `p0`: 分块起始格点索引
    /// - `p1`: 分块结束格点索引
    /// - `atm`: 该分块所属的原子索引
    /// - `rmin`: 分块内离原子核的最近距离
    /// - `rmax`: 分块内离原子核的最远距离
    #[new]
    pub fn new(p0: usize, p1: usize, atm: usize, rmin: f64, rmax: f64) -> Self {
        Self {
            core: rswfn::qmcal::ks::GridBlock {
                p0,
                p1,
                atm,
                rmin,
                rmax,
            },
        }
    }

    /// 分块起始格点索引
    #[getter]
    pub fn p0(&self) -> usize {
        self.core.p0
    }

    /// 分块结束格点索引
    #[getter]
    pub fn p1(&self) -> usize {
        self.core.p1
    }

    /// 该分块所属的原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }

    /// 分块内离原子核的最近距离
    #[getter]
    pub fn rmin(&self) -> f64 {
        self.core.rmin
    }

    /// 分块内离原子核的最远距离
    #[getter]
    pub fn rmax(&self) -> f64 {
        self.core.rmax
    }
}

/// 自旋极化电子密度
#[pyclass]
pub struct SpinRho {
    core: rswfn::qmcal::ks::SpinRho,
}

#[pymethods]
impl SpinRho {
    /// α 电子密度
    #[getter]
    pub fn rho_a(&self) -> Vec<f64> {
        self.core.rho_a.clone()
    }

    /// β 电子密度
    #[getter]
    pub fn rho_b(&self) -> Vec<f64> {
        self.core.rho_b.clone()
    }

    /// α 密度梯度 `(ngrid, 3)`
    #[getter]
    pub fn grad_a(&self, py: Python) -> Py<PyArray2<f64>> {
        grids2array(py, &self.core.grad_a)
    }

    /// β 密度梯度 `(ngrid, 3)`
    #[getter]
    pub fn grad_b(&self, py: Python) -> Py<PyArray2<f64>> {
        grids2array(py, &self.core.grad_b)
    }

    /// α 密度梯度内积 `∇ρ_α · ∇ρ_α`
    #[getter]
    pub fn saa(&self) -> Vec<f64> {
        self.core.saa.clone()
    }

    /// β 密度梯度内积 `∇ρ_β · ∇ρ_β`
    #[getter]
    pub fn sbb(&self) -> Vec<f64> {
        self.core.sbb.clone()
    }

    /// α/β 密度梯度交叉内积 `∇ρ_α · ∇ρ_β`
    #[getter]
    pub fn sab(&self) -> Vec<f64> {
        self.core.sab.clone()
    }
}

/// 计算闭壳层 Kohn-Sham 电子能量
///
/// ## 参数
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `hmat`: 核心哈密顿矩阵（NumPy 二维数组）
/// - `jmat`: 库仑矩阵（NumPy 二维数组）
/// - `kmat`: 交换矩阵（NumPy 二维数组）
/// - `exact_x`: 精确交换的混合系数
/// - `xc_eng`: 交换相关能
///
/// ## 返回
/// 电子能量（Hartree）
#[pyfunction]
pub fn rks_ele_eng(
    dmat: Bound<'_, PyArray2<f64>>,
    hmat: Bound<'_, PyArray2<f64>>,
    jmat: Bound<'_, PyArray2<f64>>,
    kmat: Bound<'_, PyArray2<f64>>,
    exact_x: f64,
    xc_eng: f64,
) -> f64 {
    let dmat = dmat.readonly().as_array().to_owned();
    let hmat = hmat.readonly().as_array().to_owned();
    let jmat = jmat.readonly().as_array().to_owned();
    let kmat = kmat.readonly().as_array().to_owned();
    rswfn::qmcal::ks::rks_ele_eng(&dmat, &hmat, &jmat, &kmat, exact_x, xc_eng)
}

/// 计算自旋极化 Kohn-Sham 电子能量
///
/// ## 参数
/// - `dmat`: 总密度矩阵（NumPy 二维数组）
/// - `hmat`: 核心哈密顿矩阵（NumPy 二维数组）
/// - `jmat`: 库仑矩阵（NumPy 二维数组）
/// - `dmat_a`: α 密度矩阵（NumPy 二维数组）
/// - `dmat_b`: β 密度矩阵（NumPy 二维数组）
/// - `kmat_a`: α 交换矩阵（NumPy 二维数组）
/// - `kmat_b`: β 交换矩阵（NumPy 二维数组）
/// - `exact_x`: 精确交换的混合系数
/// - `xc_eng`: 交换相关能
///
/// ## 返回
/// 电子能量（Hartree）
#[pyfunction]
pub fn uks_ele_eng(
    dmat: Bound<'_, PyArray2<f64>>,
    hmat: Bound<'_, PyArray2<f64>>,
    jmat: Bound<'_, PyArray2<f64>>,
    dmat_a: Bound<'_, PyArray2<f64>>,
    dmat_b: Bound<'_, PyArray2<f64>>,
    kmat_a: Bound<'_, PyArray2<f64>>,
    kmat_b: Bound<'_, PyArray2<f64>>,
    exact_x: f64,
    xc_eng: f64,
) -> f64 {
    let dmat = dmat.readonly().as_array().to_owned();
    let hmat = hmat.readonly().as_array().to_owned();
    let jmat = jmat.readonly().as_array().to_owned();
    let dmat_a = dmat_a.readonly().as_array().to_owned();
    let dmat_b = dmat_b.readonly().as_array().to_owned();
    let kmat_a = kmat_a.readonly().as_array().to_owned();
    let kmat_b = kmat_b.readonly().as_array().to_owned();
    rswfn::qmcal::ks::uks_ele_eng(
        &dmat, &hmat, &jmat, &dmat_a, &dmat_b, &kmat_a, &kmat_b, exact_x, xc_eng,
    )
}

/// 在父模块中注册 `ks` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "ks")?;
    m.add_class::<RKS>()?;
    m.add_class::<RKSResult>()?;
    m.add_class::<UKS>()?;
    m.add_class::<UKSResult>()?;
    m.add_class::<GridBlock>()?;
    m.add_class::<SpinRho>()?;
    m.add("B3LYP_EXACT_X", rswfn::qmcal::ks::B3LYP_EXACT_X)?;
    m.add_function(wrap_pyfunction!(rks_ele_eng, &m)?)?;
    m.add_function(wrap_pyfunction!(uks_ele_eng, &m)?)?;
    parent_module.add_submodule(&m)
}
