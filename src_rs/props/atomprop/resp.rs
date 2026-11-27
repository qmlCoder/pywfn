#![allow(non_snake_case)]

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

/// 静电势的类型
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Clone, PartialEq)]
pub enum EspType {
    /// 总静电势（核 + 电子）
    Total,
    /// 仅电子贡献
    Electronic,
}

impl EspType {
    fn to_core(&self) -> rswfn::props::atomprop::resp::EspType {
        match self {
            EspType::Total => rswfn::props::atomprop::resp::EspType::Total,
            EspType::Electronic => rswfn::props::atomprop::resp::EspType::Electronic,
        }
    }
}

/// RESP 拟合的参数设置
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct RespOptions {
    pub layer_scales: Vec<f64>,
    pub pt_density: f64,
    pub hyper_a1: f64,
    pub hyper_a2: f64,
    pub hyper_b: f64,
    pub max_iter: usize,
    pub conv: f64,
    pub two_stage: bool,
    pub esp_type: EspType,
}

impl Default for RespOptions {
    fn default() -> Self {
        Self {
            layer_scales: vec![],
            pt_density: 1.0,
            hyper_a1: 0.0,
            hyper_a2: 0.0,
            hyper_b: 0.0,
            max_iter: 500,
            conv: 1.0e-6,
            two_stage: true,
            esp_type: EspType::Total,
        }
    }
}

impl RespOptions {
    fn to_core(&self) -> rswfn::props::atomprop::resp::RespOptions {
        rswfn::props::atomprop::resp::RespOptions {
            layer_scales: self.layer_scales.clone(),
            pt_density: self.pt_density,
            hyper_a1: self.hyper_a1,
            hyper_a2: self.hyper_a2,
            hyper_b: self.hyper_b,
            max_iter: self.max_iter,
            conv: self.conv,
            two_stage: self.two_stage,
            esp_type: self.esp_type.to_core(),
        }
    }
}

#[pymethods]
impl RespOptions {
    /// 创建 RESP 拟合参数
    ///
    /// ## 参数
    /// - `layer_scales`: 各层格点密度的缩放系数
    /// - `pt_density`: 格点密度
    /// - `hyper_a1`: 第一层双曲惩罚系数
    /// - `hyper_a2`: 第二层双曲惩罚系数
    /// - `hyper_b`: 双曲惩罚的偏移量
    /// - `max_iter`: 最大迭代次数
    /// - `conv`: 收敛阈值
    /// - `two_stage`: 是否使用两阶段拟合
    /// - `esp_type`: 静电势类型
    #[new]
    #[pyo3(signature = (layer_scales=None, pt_density=1.0, hyper_a1=0.0, hyper_a2=0.0, hyper_b=0.0, max_iter=500, conv=1.0e-6, two_stage=true, esp_type=None))]
    pub fn new(
        layer_scales: Option<Vec<f64>>,
        pt_density: f64,
        hyper_a1: f64,
        hyper_a2: f64,
        hyper_b: f64,
        max_iter: usize,
        conv: f64,
        two_stage: bool,
        esp_type: Option<EspType>,
    ) -> Self {
        Self {
            layer_scales: layer_scales.unwrap_or_default(),
            pt_density,
            hyper_a1,
            hyper_a2,
            hyper_b,
            max_iter,
            conv,
            two_stage,
            esp_type: esp_type.unwrap_or(EspType::Total),
        }
    }

    /// 使用 rswfn 的默认参数
    #[staticmethod]
    pub fn default() -> Self {
        let core = rswfn::props::atomprop::resp::RespOptions::default();
        Self {
            layer_scales: core.layer_scales,
            pt_density: core.pt_density,
            hyper_a1: core.hyper_a1,
            hyper_a2: core.hyper_a2,
            hyper_b: core.hyper_b,
            max_iter: core.max_iter,
            conv: core.conv,
            two_stage: core.two_stage,
            esp_type: match core.esp_type {
                rswfn::props::atomprop::resp::EspType::Total => EspType::Total,
                rswfn::props::atomprop::resp::EspType::Electronic => EspType::Electronic,
            },
        }
    }

    /// 各层格点密度的缩放系数
    #[getter]
    pub fn layer_scales(&self) -> Vec<f64> {
        self.layer_scales.clone()
    }

    #[setter]
    pub fn set_layer_scales(&mut self, val: Vec<f64>) {
        self.layer_scales = val;
    }

    /// 格点密度
    #[getter]
    pub fn pt_density(&self) -> f64 {
        self.pt_density
    }

    #[setter]
    pub fn set_pt_density(&mut self, val: f64) {
        self.pt_density = val;
    }

    /// 第一层双曲惩罚系数
    #[getter]
    pub fn hyper_a1(&self) -> f64 {
        self.hyper_a1
    }

    #[setter]
    pub fn set_hyper_a1(&mut self, val: f64) {
        self.hyper_a1 = val;
    }

    /// 第二层双曲惩罚系数
    #[getter]
    pub fn hyper_a2(&self) -> f64 {
        self.hyper_a2
    }

    #[setter]
    pub fn set_hyper_a2(&mut self, val: f64) {
        self.hyper_a2 = val;
    }

    /// 双曲惩罚的偏移量
    #[getter]
    pub fn hyper_b(&self) -> f64 {
        self.hyper_b
    }

    #[setter]
    pub fn set_hyper_b(&mut self, val: f64) {
        self.hyper_b = val;
    }

    /// 最大迭代次数
    #[getter]
    pub fn max_iter(&self) -> usize {
        self.max_iter
    }

    #[setter]
    pub fn set_max_iter(&mut self, val: usize) {
        self.max_iter = val;
    }

    /// 收敛阈值
    #[getter]
    pub fn conv(&self) -> f64 {
        self.conv
    }

    #[setter]
    pub fn set_conv(&mut self, val: f64) {
        self.conv = val;
    }

    /// 是否使用两阶段拟合
    #[getter]
    pub fn two_stage(&self) -> bool {
        self.two_stage
    }

    #[setter]
    pub fn set_two_stage(&mut self, val: bool) {
        self.two_stage = val;
    }

    /// 静电势类型
    #[getter]
    pub fn esp_type(&self) -> EspType {
        self.esp_type.clone()
    }

    #[setter]
    pub fn set_esp_type(&mut self, val: EspType) {
        self.esp_type = val;
    }
}

/// RESP 拟合结果
#[pyclass]
#[derive(Clone)]
pub struct RespResult {
    pub charges: Vec<f64>,
    pub stage1: Vec<f64>,
    pub n_points: usize,
    pub rmse: f64,
    pub rrmse: f64,
    pub stage1_iter: usize,
    pub stage2_iter: usize,
    pub stage2: bool,
    pub converged: bool,
}

impl RespResult {
    fn from_core(core: rswfn::props::atomprop::resp::RespResult) -> Self {
        Self {
            charges: core.charges,
            stage1: core.stage1,
            n_points: core.n_points,
            rmse: core.rmse,
            rrmse: core.rrmse,
            stage1_iter: core.stage1_iter,
            stage2_iter: core.stage2_iter,
            stage2: core.stage2,
            converged: core.converged,
        }
    }
}

#[pymethods]
impl RespResult {
    /// 拟合得到的电荷
    #[getter]
    pub fn charges(&self) -> Vec<f64> {
        self.charges.clone()
    }

    /// 第一阶段得到的电荷
    #[getter]
    pub fn stage1(&self) -> Vec<f64> {
        self.stage1.clone()
    }

    /// 拟合使用的格点数量
    #[getter]
    pub fn n_points(&self) -> usize {
        self.n_points
    }

    /// 均方根误差
    #[getter]
    pub fn rmse(&self) -> f64 {
        self.rmse
    }

    /// 相对均方根误差
    #[getter]
    pub fn rrmse(&self) -> f64 {
        self.rrmse
    }

    /// 第一阶段的迭代次数
    #[getter]
    pub fn stage1_iter(&self) -> usize {
        self.stage1_iter
    }

    /// 第二阶段的迭代次数
    #[getter]
    pub fn stage2_iter(&self) -> usize {
        self.stage2_iter
    }

    /// 是否进行了第二阶段拟合
    #[getter]
    pub fn stage2(&self) -> bool {
        self.stage2
    }

    /// 是否收敛
    #[getter]
    pub fn converged(&self) -> bool {
        self.converged
    }
}

/// RESP 电荷计算器
#[pyclass]
pub struct Calculator {
    pub mole: Mole,
    pub atms: Vec<usize>,
}

impl Calculator {
    /// 创建底层的 RESP 计算器
    fn caler(&self) -> rswfn::props::atomprop::resp::Calculator<'_> {
        let mut caler = rswfn::props::atomprop::resp::Calculator::new(&self.mole.core);
        caler.atms = &self.atms;
        caler
    }
}

#[pymethods]
impl Calculator {
    /// 创建一个新的 RESP 计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `atms`: 可选的需要拟合的原子索引列表，默认全部原子
    #[new]
    #[pyo3(signature = (mole, atms=None))]
    pub fn new(mole: Mole, atms: Option<Vec<usize>>) -> Self {
        let atms = atms.unwrap_or_else(|| mole.core.atoms().idxs().clone());
        Self { mole, atms }
    }

    /// 使用默认参数计算 RESP 电荷
    ///
    /// ## 返回
    /// RESP 拟合结果
    pub fn resp(&self) -> PyResult<RespResult> {
        self.caler()
            .resp()
            .map(RespResult::from_core)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))
    }

    /// 使用指定参数计算 RESP 电荷
    ///
    /// ## 参数
    /// - `opts`: 拟合参数
    ///
    /// ## 返回
    /// RESP 拟合结果
    pub fn resp_with(&self, opts: RespOptions) -> PyResult<RespResult> {
        self.caler()
            .resp_with(&opts.to_core())
            .map(RespResult::from_core)
            .map_err(|e| PyRuntimeError::new_err(format!("{e}")))
    }

    /// 生成用于拟合的格点
    ///
    /// ## 参数
    /// - `opts`: 拟合参数
    ///
    /// ## 返回
    /// 格点坐标列表
    pub fn mk_points(&self, opts: RespOptions) -> Vec<[f64; 3]> {
        self.caler().mk_points(&opts.to_core())
    }

    /// 计算指定格点上的静电势
    ///
    /// ## 参数
    /// - `pts`: 格点坐标列表
    /// - `esp_type`: 静电势类型
    ///
    /// ## 返回
    /// 每个格点上的静电势列表
    pub fn esp_values(&self, pts: Vec<[f64; 3]>, esp_type: EspType) -> Vec<f64> {
        self.caler().esp_values(&pts, esp_type.to_core())
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "resp")?;
    m.add_class::<EspType>()?;
    m.add_class::<RespOptions>()?;
    m.add_class::<RespResult>()?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}
