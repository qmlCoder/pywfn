use crate::base::Mole;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

/// 返回每个原子轨道所属的原子索引
///
/// ## 参数
/// - `mole`: 分子对象
/// - `form`: 基函数形式（"sph" 或 "car"）
///
/// ## 返回
/// 长度为原子轨道数的列表，每个元素是该轨道所属的原子索引
#[pyfunction]
pub fn ao_owners(mole: Mole, form: String) -> Vec<usize> {
    rswfn::qmcal::guess::ao_owners(&mole.core, &form)
}

/// 是否启用 SAD（叠加原子密度）初猜
///
/// ## 返回
/// 启用时返回 `True`，否则返回 `False`
#[pyfunction]
pub fn use_sad_enabled() -> bool {
    rswfn::qmcal::guess::use_sad_enabled()
}

/// 用叠加原子密度（SAD）方法构造初始密度矩阵
///
/// ## 参数
/// - `smat`: 重叠矩阵（NumPy 二维数组）
/// - `hmat`: 核心哈密顿矩阵（NumPy 二维数组）
/// - `ao_atm`: 每个原子轨道所属的原子索引
/// - `nucs`: 各原子的原子序数
/// - `nele`: 电子数
///
/// ## 返回
/// 初始密度矩阵（NumPy 二维数组）；无法构造时返回 `None`
#[pyfunction]
pub fn sad_dmat(
    py: Python,
    smat: Bound<'_, PyArray2<f64>>,
    hmat: Bound<'_, PyArray2<f64>>,
    ao_atm: Vec<usize>,
    nucs: Vec<usize>,
    nele: usize,
) -> Option<Py<PyArray2<f64>>> {
    let smat = smat.readonly().as_array().to_owned();
    let hmat = hmat.readonly().as_array().to_owned();
    rswfn::qmcal::guess::sad_dmat(&smat, &hmat, &ao_atm, &nucs, nele)
        .map(|dmat| dmat.into_pyarray(py).unbind())
}

/// 在父模块中注册 `guess` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "guess")?;
    m.add_function(wrap_pyfunction!(ao_owners, &m)?)?;
    m.add_function(wrap_pyfunction!(use_sad_enabled, &m)?)?;
    m.add_function(wrap_pyfunction!(sad_dmat, &m)?)?;
    parent_module.add_submodule(&m)
}
