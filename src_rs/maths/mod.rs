use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

pub mod march;
pub mod matrix;

/// 拉格朗日插值
///
/// 给定数据点 (xs, ys)，在目标点 ts 处进行拉格朗日插值。
///
/// ## 参数
/// - `xs`: 已知点的 x 坐标列表
/// - `ys`: 已知点的 y 坐标列表
/// - `ts`: 需要插值的目标点列表
///
/// ## 返回
/// 返回插值结果列表
#[pyfunction]
pub fn lag_intpol(xs: Vec<f64>, ys: Vec<f64>, ts: Vec<f64>) -> Vec<f64> {
    ts.iter()
        .map(|t| rswfn::maths::lag_intpol(&xs, &ys, *t))
        .collect()
}

/// 计算矩阵的平方根逆矩阵 P^(-1/2)
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `mat`: 输入的对称正定矩阵（二维 NumPy 数组）
///
/// ## 返回
/// 返回矩阵的 P^(-1/2)（二维 NumPy 数组）
#[pyfunction]
pub fn matrix_phalf(py: Python, mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let mat = mat.readonly().as_array().to_owned();
    let phalf = rswfn::maths::matrix_phalf(&mat);
    phalf.into_pyarray(py).unbind()
}

/// 计算矩阵的负平方根逆矩阵 P^(-1/2)（另一种约定）
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `mat`: 输入的对称正定矩阵（二维 NumPy 数组）
///
/// ## 返回
/// 返回矩阵的 N^(-1/2)（二维 NumPy 数组）
#[pyfunction]
pub fn matrix_nhalf(py: Python, mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let mat = mat.readonly().as_array().to_owned();
    let nhalf = rswfn::maths::matrix_nhalf(&mat);
    nhalf.into_pyarray(py).unbind()
}

/// 从分子轨道系数和占据数计算密度矩阵
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `cmat`: 分子轨道系数矩阵（二维 NumPy 数组）
/// - `obt_occs`: 分子轨道占据数列表
///
/// ## 返回
/// 返回密度矩阵（二维 NumPy 数组）
#[pyfunction]
pub fn calc_dmat(
    py: Python,
    cmat: Bound<'_, PyArray2<f64>>,
    obt_occs: Vec<f64>,
) -> Py<PyArray2<f64>> {
    let cmat = cmat.readonly().as_array().to_owned();
    let dmat = rswfn::maths::calc_dmat(&cmat, &obt_occs);
    dmat.into_pyarray(py).unbind()
}

/// 计算阶乘 n!
///
/// ## 参数
/// - `n`: 非负整数
///
/// ## 返回
/// 返回 n! 的值
#[pyfunction]
pub fn fac(n: isize) -> isize {
    rswfn::maths::fac(n)
}

/// 计算双阶乘 n!!
///
/// ## 参数
/// - `n`: 整数
///
/// ## 返回
/// 返回 n!! 的值
#[pyfunction]
pub fn fac2(n: isize) -> isize {
    rswfn::maths::fac2(n)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "maths")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(lag_intpol, &m)?)?; //
    m.add_function(wrap_pyfunction!(calc_dmat, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_phalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_nhalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac2, &m)?)?; //
    march::register_module(&m)?;
    matrix::register_module(&m)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
