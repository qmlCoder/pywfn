use crate::base::{Basis, Coefs, Mole};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

/// RHF（闭壳层 Hartree-Fock）计算结果
#[pyclass]
pub struct RHFResult {
    core: rswfn::qmcal::hf::RHFResult,
}

#[pymethods]
impl RHFResult {
    /// 电子能量（Hartree）
    #[getter]
    pub fn ele_eng(&self) -> f64 {
        self.core.ele_eng
    }

    /// 核-核排斥能（Hartree）
    #[getter]
    pub fn nuc_eng(&self) -> f64 {
        self.core.nuc_eng
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

    /// 分子总能量（电子能量 + 核排斥能）
    pub fn total_eng(&self) -> f64 {
        self.core.total_eng()
    }
}

/// RHF（闭壳层 Hartree-Fock）计算器
///
/// 提供闭壳层自洽场（SCF）计算功能，用于求解分子的 Hartree-Fock 方程。
#[pyclass]
pub struct RHF {
    mole: Mole,
    opt_int_form: Option<String>,
    opt_max_iter: Option<usize>,
    opt_tol: Option<f64>,
    opt_screen: Option<f64>,
    opt_max_mem: Option<usize>,
}

impl RHF {
    fn caler(&self) -> rswfn::qmcal::hf::RHF<'_> {
        let mut caler = rswfn::qmcal::hf::RHF::new(&self.mole.core);
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
        caler
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
        let rhf = Self {
            mole,
            opt_int_form: None,
            opt_max_iter: None,
            opt_tol: None,
            opt_screen: None,
            opt_max_mem: None,
        };
        Ok(rhf)
    }

    /// 执行 RHF 自洽场计算。
    ///
    /// 迭代求解 Hartree-Fock 方程，获得收敛的分子轨道和能量。
    /// 双电子积分由 rswfn 内部的壳层块引擎自动计算（直接 SCF）。
    ///
    /// ## 返回
    /// RHF 计算结果对象
    pub fn run(&self) -> RHFResult {
        RHFResult {
            core: self.caler().run(),
        }
    }

    /// 核-核排斥能（Hartree）
    pub fn nuc_eng(&self) -> f64 {
        self.caler().nuc_eng()
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

    /// 最大 SCF 迭代次数
    #[getter]
    pub fn max_iter(&self) -> usize {
        self.caler().max_iter
    }
    #[setter]
    pub fn set_max_iter(&mut self, value: usize) {
        self.opt_max_iter = Some(value);
    }

    /// SCF 收敛阈值（能量变化）
    #[getter]
    pub fn tol(&self) -> f64 {
        self.caler().tol
    }
    #[setter]
    pub fn set_tol(&mut self, value: f64) {
        self.opt_tol = Some(value);
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
}

/// UHF（无限制 Hartree-Fock）计算器
#[pyclass]
pub struct UHF {
    mole: Mole,
}

impl UHF {
    fn caler(&self) -> rswfn::qmcal::hf::UHF<'_> {
        rswfn::qmcal::hf::UHF::new(&self.mole.core)
    }
}

#[pymethods]
impl UHF {
    /// 创建 UHF 计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 由 α/β 轨道系数构造密度矩阵
    ///
    /// ## 参数
    /// - `nele`: 二元组，α 与 β 电子数 `[nela, nelb]`
    /// - `cmat_a`: α 轨道系数矩阵（NumPy 二维数组）
    /// - `cmat_b`: β 轨道系数矩阵（NumPy 二维数组）
    ///
    /// ## 返回
    /// 二元组 `(α 密度矩阵, β 密度矩阵)`
    pub fn calc_udmat(
        &self,
        py: Python,
        nele: [usize; 2],
        cmat_a: Bound<'_, PyArray2<f64>>,
        cmat_b: Bound<'_, PyArray2<f64>>,
    ) -> (Py<PyArray2<f64>>, Py<PyArray2<f64>>) {
        let coe_a = cmat_a.readonly().as_array().to_owned();
        let coe_b = cmat_b.readonly().as_array().to_owned();
        let cmats = [coe_a, coe_b];
        let [dmat_a, dmat_b] = self.caler().calc_udmat(nele, &cmats);
        (
            dmat_a.into_pyarray(py).unbind(),
            dmat_b.into_pyarray(py).unbind(),
        )
    }

    /// 执行无限制 Hartree-Fock 计算
    ///
    /// ## 返回
    /// 三元组 `(α/β 电子能量, (α 轨道能量, β 轨道能量), (α 轨道系数, β 轨道系数))`
    pub fn calc_uhf(
        &self,
        py: Python,
    ) -> (
        [f64; 2],
        (Py<PyArray1<f64>>, Py<PyArray1<f64>>),
        (Py<PyArray2<f64>>, Py<PyArray2<f64>>),
    ) {
        let (engs, obts, cmats) = self.caler().calc_uhf();
        let [eng_a, eng_b] = engs;
        let [obt_a, obt_b] = obts;
        let [coef_a, coef_b] = cmats;
        (
            [eng_a, eng_b],
            (
                obt_a.into_pyarray(py).unbind(),
                obt_b.into_pyarray(py).unbind(),
            ),
            (
                coef_a.into_pyarray(py).unbind(),
                coef_b.into_pyarray(py).unbind(),
            ),
        )
    }
}

/// 计算闭壳层电子能量
///
/// ## 参数
/// - `dmat`: 密度矩阵（NumPy 二维数组）
/// - `hmat`: 核心哈密顿矩阵（NumPy 二维数组）
/// - `fmat`: Fock 矩阵（NumPy 二维数组）
///
/// ## 返回
/// 电子能量（Hartree）
#[pyfunction]
pub fn rhf_ele_eng(
    dmat: Bound<'_, PyArray2<f64>>,
    hmat: Bound<'_, PyArray2<f64>>,
    fmat: Bound<'_, PyArray2<f64>>,
) -> f64 {
    let dmat = dmat.readonly().as_array().to_owned();
    let hmat = hmat.readonly().as_array().to_owned();
    let fmat = fmat.readonly().as_array().to_owned();
    rswfn::qmcal::hf::rhf_ele_eng(&dmat, &hmat, &fmat)
}

/// 计算开壳层（无限制）电子能量
///
/// ## 参数
/// - `dmat_a`: α 密度矩阵（NumPy 二维数组）
/// - `dmat_b`: β 密度矩阵（NumPy 二维数组）
/// - `tmat`: 动能矩阵（NumPy 二维数组）
/// - `vmat`: 核吸引矩阵（NumPy 二维数组）
/// - `fmat_a`: α Fock 矩阵（NumPy 二维数组）
/// - `fmat_b`: β Fock 矩阵（NumPy 二维数组）
///
/// ## 返回
/// 二元组 `[α 电子能量, β 电子能量]`
#[pyfunction]
pub fn uhf_ele_eng(
    dmat_a: Bound<'_, PyArray2<f64>>,
    dmat_b: Bound<'_, PyArray2<f64>>,
    tmat: Bound<'_, PyArray2<f64>>,
    vmat: Bound<'_, PyArray2<f64>>,
    fmat_a: Bound<'_, PyArray2<f64>>,
    fmat_b: Bound<'_, PyArray2<f64>>,
) -> [f64; 2] {
    let dmat_a = dmat_a.readonly().as_array().to_owned();
    let dmat_b = dmat_b.readonly().as_array().to_owned();
    let tmat = tmat.readonly().as_array().to_owned();
    let vmat = vmat.readonly().as_array().to_owned();
    let fmat_a = fmat_a.readonly().as_array().to_owned();
    let fmat_b = fmat_b.readonly().as_array().to_owned();
    let dmat = [dmat_a, dmat_b];
    let fmat = [fmat_a, fmat_b];
    rswfn::qmcal::hf::uhf_ele_eng(&dmat, &tmat, &vmat, &fmat)
}

/// 在父模块中注册 `hf` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "hf")?;
    m.add_class::<RHF>()?;
    m.add_class::<RHFResult>()?;
    m.add_class::<UHF>()?;
    m.add_function(wrap_pyfunction!(rhf_ele_eng, &m)?)?;
    m.add_function(wrap_pyfunction!(uhf_ele_eng, &m)?)?;
    parent_module.add_submodule(&m)
}
