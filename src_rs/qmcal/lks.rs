use std::collections::HashMap;

use crate::base::{Mole, Stm};
use crate::qmcal::ks::grids2array;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

/// 局部 π 电子 Kohn-Sham（LKS）计算结果
#[pyclass]
pub struct LKSResult {
    core: rswfn::qmcal::lks::LKSResult,
}

#[pymethods]
impl LKSResult {
    /// 电子能量（含交换相关能），单位 Hartree
    #[getter]
    pub fn ele_eng(&self) -> f64 {
        self.core.ele_eng
    }
    /// 交换相关能（Hartree）
    #[getter]
    pub fn xc_eng(&self) -> f64 {
        self.core.xc_eng
    }
    /// 核-核排斥能（Hartree）
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
    /// π 电子数
    #[getter]
    pub fn npi_elec(&self) -> usize {
        self.core.npi_elec
    }
    /// π 轨道数
    #[getter]
    pub fn npi(&self) -> usize {
        self.core.npi
    }
    /// 原子轨道数
    #[getter]
    pub fn nato(&self) -> usize {
        self.core.nato
    }
    /// 分子轨道能量
    #[getter]
    pub fn obt_engs(&self) -> Vec<f64> {
        self.core.obt_engs.clone()
    }
    /// π 原子轨道索引
    #[getter]
    pub fn pi_atos(&self) -> Vec<usize> {
        self.core.pi_atos.clone()
    }
    /// 各 π 原子轨道对应的局域坐标系
    #[getter]
    pub fn stms(&self) -> HashMap<usize, Stm> {
        self.core
            .stms
            .iter()
            .map(|(atm, stm)| (*atm, Stm { core: stm.clone() }))
            .collect()
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
    /// π 轨道系数矩阵（nato × npi）
    pub fn cmat(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.cmat.clone().into_pyarray(py).unbind()
    }
    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.core.total_eng()
    }
    /// 占据轨道数
    pub fn nocc(&self) -> usize {
        self.core.nocc()
    }
    /// 把局域 π 轨道系数展开为全基组的系数矩阵
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 展开后的系数矩阵（NumPy 二维数组）
    pub fn expand_coefs(&self, py: Python, mole: Mole) -> Py<PyArray2<f64>> {
        self.core
            .expand_coefs(&mole.core)
            .into_pyarray(py)
            .unbind()
    }
    /// 返回各局域 π 轨道的标签
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 标签字符串列表
    pub fn local_labels(&self, mole: Mole) -> Vec<String> {
        self.core.local_labels(&mole.core)
    }
    /// 返回局域 π 轨道的表格化文本
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 轨道信息表格字符串
    pub fn mo_table_local(&self, mole: Mole) -> String {
        self.core.mo_table_local(&mole.core)
    }
    /// 把局域 π 计算结果还原成分子对象（写入轨道系数）
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 带有本次计算轨道系数的新分子对象
    pub fn to_mole(&self, mole: Mole) -> Mole {
        Mole {
            core: self.core.to_mole(&mole.core),
        }
    }
}

/// 局部 π 电子 Kohn-Sham（LKS）计算器
#[pyclass]
pub struct LKS {
    pub mole: Mole,
    pub func: String,
    pub grid: (usize, usize),
    pub npi_elec: Option<usize>,
    pub opt_max_iter: Option<usize>,
    pub opt_tol: Option<f64>,
    pub opt_screen: Option<f64>,
    pub opt_max_mem: Option<usize>,
}

impl LKS {
    fn build(
        mole: Mole,
        func: String,
        grid: (usize, usize),
        npi_elec: Option<usize>,
    ) -> Self {
        Self {
            mole,
            func,
            grid,
            npi_elec,
            opt_max_iter: None,
            opt_tol: None,
            opt_screen: None,
            opt_max_mem: None,
        }
    }

    fn caler(&self) -> rswfn::qmcal::lks::LKS<'_> {
        let func = rswfn::rsxcs::XcFunc::from_name(&self.func).expect("未知泛函");
        let mut caler = rswfn::qmcal::lks::LKS::with_grid(
            &self.mole.core,
            self.grid.0,
            self.grid.1,
            func,
            self.npi_elec,
        );
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
impl LKS {
    /// 按泛函名创建 LKS 计算器（默认格点 100×74）。
    #[new]
    pub fn new(mole: Mole, func: String) -> PyResult<Self> {
        Self::check(&func)?;
        Ok(Self::build(mole, func, (100, 74), None))
    }

    /// 按泛函名创建 LKS 计算器（默认格点 100×74）；泛函名无效时返回 `None`。
    #[staticmethod]
    pub fn from_name(mole: Mole, name: String) -> Option<Self> {
        rswfn::rsxcs::XcFunc::from_name(&name)?;
        Some(Self::build(mole, name, (100, 74), None))
    }

    /// 用指定格点与 π 电子数创建计算器。
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `radial`: 径向格点数
    /// - `angular`: 球面格点数
    /// - `func`: 泛函名（如 "b3lyp"）
    /// - `npi_elec`: π 电子数；为 `None` 时由分子自动推断
    #[staticmethod]
    pub fn with_grid(
        mole: Mole,
        radial: usize,
        angular: usize,
        func: String,
        npi_elec: Option<usize>,
    ) -> PyResult<Self> {
        Self::check(&func)?;
        Ok(Self::build(mole, func, (radial, angular), npi_elec))
    }

    /// 执行自洽场计算。
    pub fn run(&self) -> LKSResult {
        LKSResult {
            core: self.caler().run(),
        }
    }

    /// 核-核排斥能（Hartree）
    pub fn nuc_eng(&self) -> f64 {
        self.caler().nuc_eng()
    }

    /// 泛函名称
    #[getter]
    pub fn func(&self) -> String {
        self.func.clone()
    }

    /// 径向格点数
    #[getter]
    pub fn radial(&self) -> usize {
        self.grid.0
    }

    /// 球面格点数
    #[getter]
    pub fn angular(&self) -> usize {
        self.grid.1
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

    /// π 电子数；`None` 表示由分子自动推断
    #[getter]
    pub fn npi_elec(&self) -> Option<usize> {
        self.npi_elec
    }
    #[setter]
    pub fn set_npi_elec(&mut self, value: Option<usize>) {
        self.npi_elec = value;
    }

    /// 各 π 原子轨道对应的局域坐标系
    #[getter]
    pub fn stms(&self) -> HashMap<usize, Stm> {
        let caler = self.caler();
        caler
            .stms
            .clone()
            .into_iter()
            .map(|(atm, stm)| (atm, Stm { core: stm }))
            .collect()
    }

    /// π 原子轨道索引
    #[getter]
    pub fn pi_atos(&self) -> Vec<usize> {
        self.caler().pi_atos.clone()
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

    /// π 轨道数
    #[getter]
    pub fn npi(&self) -> usize {
        self.caler().npi
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

    /// π 基函数之间的变换矩阵列表
    #[getter]
    pub fn ts(&self, py: Python) -> Vec<Py<PyArray2<f64>>> {
        let caler = self.caler();
        caler
            .ts
            .clone()
            .into_iter()
            .map(|mat| mat.into_pyarray(py).unbind())
            .collect()
    }

    /// 各 π 片段包含的原子轨道索引对
    #[getter]
    pub fn pi_comp(&self) -> Vec<Vec<(usize, usize)>> {
        self.caler().pi_comp.clone()
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
}

/// 在父模块中注册 `lks` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "lks")?;
    m.add_class::<LKS>()?;
    m.add_class::<LKSResult>()?;
    parent_module.add_submodule(&m)
}
