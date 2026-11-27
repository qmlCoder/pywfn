use crate::base::Mole;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

/// dHMO 计算结果
#[pyclass]
pub struct DmoRes {
    core: rswfn::qmcal::dhmo::DmoRes,
}

#[pymethods]
impl DmoRes {
    /// 原子轨道列表，每个元素为 (原子索引, 壳层索引, 轨道符号, 方向向量 [x,y,z])
    #[getter]
    pub fn atos(&self) -> Vec<(usize, usize, String, Vec<f64>)> {
        self.core
            .atos
            .iter()
            .map(|(atm, shl, sym, dir)| (*atm, *shl, sym.clone(), dir.to_vec()))
            .collect()
    }
    /// 键连矩阵
    pub fn bm(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.BM.clone().into_pyarray(py).unbind()
    }
    /// 轨道系数矩阵
    pub fn cm(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.CM.clone().into_pyarray(py).unbind()
    }
    /// 键连矩阵（对应上游字段 `BM`）
    #[allow(non_snake_case)]
    #[getter]
    pub fn BM(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.BM.clone().into_pyarray(py).unbind()
    }
    /// 轨道系数矩阵（对应上游字段 `CM`）
    #[allow(non_snake_case)]
    #[getter]
    pub fn CM(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.CM.clone().into_pyarray(py).unbind()
    }
    /// 分子轨道能量
    #[getter]
    pub fn es(&self) -> Vec<f64> {
        self.core.es.clone()
    }
}

/// 计算 sp2 杂化原子的有效核电荷数
///
/// ## 参数
/// - `mole`: 分子对象
/// - `atm`: 原子索引
///
/// ## 返回
/// 该原子的有效核电荷数
#[pyfunction]
pub fn effective_charge(mole: Mole, atm: usize) -> f64 {
    rswfn::qmcal::dhmo::effective_charge(&mole.core, atm)
}

/// 计算两个原子间主量子数最高的 Slater 型 p 轨道的重叠积分
///
/// ## 参数
/// - `mole`: 分子对象
/// - `ai`: 第一个原子索引
/// - `aj`: 第二个原子索引
///
/// ## 返回
/// 重叠积分值
#[pyfunction]
pub fn overlap_integral(mole: Mole, ai: usize, aj: usize) -> f64 {
    rswfn::qmcal::dhmo::overlap_integral(&mole.core, ai, aj)
}

/// 计算原子的自由价
///
/// ## 参数
/// - `mole`: 分子对象
/// - `atm`: 原子索引
/// - `sum_order`: 该原子所有键级之和
///
/// ## 返回
/// 自由价
#[pyfunction]
pub fn free_valence(mole: Mole, atm: usize, sum_order: f64) -> f64 {
    rswfn::qmcal::dhmo::free_valence(&mole.core, atm, sum_order)
}

/// 执行 dHMO 计算，返回分子结构与结果
///
/// ## 参数
/// - `path`: 输入文件（.dmi）路径
///
/// ## 返回
/// 二元组 `(分子对象, dHMO 计算结果)`
#[pyfunction]
pub fn dhmo_calc(path: String) -> (Mole, DmoRes) {
    let (mole, res) = rswfn::qmcal::dhmo::dhmo_calc(path);
    (Mole { core: *mole }, DmoRes { core: res })
}

/// 计算并输出 dHMO 结果到 .dmo 文件
///
/// ## 参数
/// - `dmi_path`: 输入文件（.dmi）路径
/// - `dmo_path`: 输出文件（.dmo）路径；为 `None` 时按默认规则命名
#[pyfunction]
pub fn output(dmi_path: &str, dmo_path: Option<String>) {
    rswfn::qmcal::dhmo::output(dmi_path, dmo_path)
}

/// 由 Hückel 哈密顿矩阵求 dHMO 的轨道能量与相关量
///
/// ## 参数
/// - `mole`: 分子对象
/// - `bm`: 键连矩阵（NumPy 二维数组）
/// - `cm`: 轨道系数矩阵（NumPy 二维数组）
/// - `ato_atms`: 每个原子轨道所属的原子索引
/// - `nele`: 电子数
///
/// ## 返回
/// 三元组 `(密度矩阵, 轨道能量列表, 电子能量)`
#[pyfunction]
#[allow(non_snake_case)]
pub fn HMO_Analysis(
    py: Python,
    mole: Mole,
    bm: Bound<'_, PyArray2<f64>>,
    cm: Bound<'_, PyArray2<f64>>,
    ato_atms: Vec<usize>,
    nele: f64,
) -> PyResult<(Py<PyArray2<f64>>, Vec<f64>, f64)> {
    let bm = bm.readonly().as_array().to_owned();
    let cm = cm.readonly().as_array().to_owned();
    let (dmat, es, eng) = rswfn::qmcal::dhmo::HMO_Analysis(&mole.core, &bm, &cm, &ato_atms, nele);
    Ok((dmat.into_pyarray(py).unbind(), es, eng))
}

/// 在父模块中注册 `dhmo` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "dhmo")?;
    m.add_class::<DmoRes>()?;
    m.add_function(wrap_pyfunction!(effective_charge, &m)?)?;
    m.add_function(wrap_pyfunction!(overlap_integral, &m)?)?;
    m.add_function(wrap_pyfunction!(free_valence, &m)?)?;
    m.add_function(wrap_pyfunction!(dhmo_calc, &m)?)?;
    m.add_function(wrap_pyfunction!(output, &m)?)?;
    m.add_function(wrap_pyfunction!(HMO_Analysis, &m)?)?;
    parent_module.add_submodule(&m)
}
