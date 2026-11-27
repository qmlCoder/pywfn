use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

use crate::base::Vec3;

pub mod march;
pub mod matrix;
pub mod optim;

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

/// 获取列表中的最大值及其索引
///
/// ## 参数
/// - `vals`: 数值列表
///
/// ## 返回
/// 返回 (索引, 最大值)
#[pyfunction]
pub fn arg_max(vals: Vec<f64>) -> (usize, f64) {
    rswfn::maths::arg_max(&vals)
}

/// 获取列表中的最小值及其索引
///
/// ## 参数
/// - `vals`: 数值列表
///
/// ## 返回
/// 返回 (索引, 最小值)
#[pyfunction]
pub fn arg_min(vals: Vec<f64>) -> (usize, f64) {
    rswfn::maths::arg_min(&vals)
}

/// 按玻尔兹曼分布计算各状态的概率
///
/// ## 参数
/// - `ens`: 各状态的能量列表
/// - `t`: 温度
///
/// ## 返回
/// 返回各状态的概率列表
#[pyfunction]
pub fn bolzman_dist(ens: Vec<f64>, t: f64) -> Vec<f64> {
    rswfn::maths::bolzman_dist(&ens, t)
}

/// 计算两点之间的距离
///
/// ## 参数
/// - `p1`: 第一个点坐标
/// - `p2`: 第二个点坐标
#[pyfunction]
pub fn calc_dist(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    rswfn::maths::calc_dist(&p1, &p2)
}

/// 计算矩阵的行列式
///
/// ## 参数
/// - `mat`: 方阵
#[pyfunction]
pub fn det(mat: Bound<'_, PyArray2<f64>>) -> f64 {
    let mat = mat.readonly().as_array().to_owned();
    rswfn::maths::det(&mat)
}

/// 对对称矩阵做特征分解
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `mat`: 对称矩阵
///
/// ## 返回
/// 返回 (特征值数组, 特征向量矩阵)
#[pyfunction]
pub fn eigh(
    py: Python,
    mat: Bound<'_, PyArray2<f64>>,
) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>) {
    let mat = mat.readonly().as_array().to_owned();
    let (vals, vecs) = rswfn::maths::eigh(&mat);
    (vals.into_pyarray(py).unbind(), vecs.into_pyarray(py).unbind())
}

/// 由三个点确定一个平面
///
/// ## 参数
/// - `p1`: 平面上的第一个点
/// - `p2`: 平面上的第二个点
/// - `p3`: 平面上的第三个点
///
/// ## 返回
/// 返回平面方程 ax + by + cz + d = 0 的系数 [a, b, c, d]
#[pyfunction]
pub fn get_plan_by_3points(p1: Vec3, p2: Vec3, p3: Vec3) -> [f64; 4] {
    rswfn::maths::get_plan_by_3points(&p1.core, &p2.core, &p3.core)
}

/// 计算矩阵的逆矩阵
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `mat`: 方阵
///
/// ## 返回
/// 返回逆矩阵
#[pyfunction]
pub fn inv(py: Python, mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
    let mat = mat.readonly().as_array().to_owned();
    rswfn::maths::inv(&mat).into_pyarray(py).unbind()
}

/// 返回与 b 同号的 a
///
/// ## 参数
/// - `a`: 数值
/// - `b`: 决定符号的数值
#[pyfunction]
pub fn sign(a: f64, b: f64) -> f64 {
    rswfn::maths::sign(a, b)
}

/// 计算两个向量的叉积
///
/// ## 参数
/// - `vec1`: 第一个向量
/// - `vec2`: 第二个向量
#[pyfunction]
pub fn vec_cross(vec1: [f64; 3], vec2: [f64; 3]) -> [f64; 3] {
    rswfn::maths::vec_cross(&vec1, &vec2)
}

/// 计算两个向量的距离
///
/// ## 参数
/// - `vec1`: 第一个向量
/// - `vec2`: 第二个向量
#[pyfunction]
pub fn vec_dist(vec1: [f64; 3], vec2: [f64; 3]) -> f64 {
    rswfn::maths::vec_dist(&vec1, &vec2)
}

/// 计算向量的长度
///
/// ## 参数
/// - `vec`: 向量
#[pyfunction]
pub fn vec_len(vec: [f64; 3]) -> f64 {
    rswfn::maths::vec_len(&vec)
}

/// 将向量归一化
///
/// ## 参数
/// - `vec`: 向量
#[pyfunction]
pub fn vec_norm(vec: [f64; 3]) -> [f64; 3] {
    rswfn::maths::vec_norm(&vec)
}

/// 计算向量在另一个向量方向上的投影
///
/// ## 参数
/// - `vec1`: 待投影的向量
/// - `vec2`: 投影方向
#[pyfunction]
pub fn vec_project(vec1: [f64; 3], vec2: [f64; 3]) -> PyResult<[f64; 3]> {
    rswfn::maths::vec_project(vec1, vec2).map_err(|e| pyo3::exceptions::PyValueError::new_err(e))
}

/// 计算向量的平方和
///
/// ## 参数
/// - `nums`: 数值列表
#[pyfunction]
pub fn vec_psd(nums: Vec<f64>) -> f64 {
    rswfn::maths::vec_psd(&nums)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "maths")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(lag_intpol, &m)?)?; //
    m.add_function(wrap_pyfunction!(calc_dmat, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_phalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(matrix_nhalf, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac, &m)?)?; //
    m.add_function(wrap_pyfunction!(fac2, &m)?)?; //
    m.add_function(wrap_pyfunction!(arg_max, &m)?)?;
    m.add_function(wrap_pyfunction!(arg_min, &m)?)?;
    m.add_function(wrap_pyfunction!(bolzman_dist, &m)?)?;
    m.add_function(wrap_pyfunction!(calc_dist, &m)?)?;
    m.add_function(wrap_pyfunction!(det, &m)?)?;
    m.add_function(wrap_pyfunction!(eigh, &m)?)?;
    m.add_function(wrap_pyfunction!(get_plan_by_3points, &m)?)?;
    m.add_function(wrap_pyfunction!(inv, &m)?)?;
    m.add_function(wrap_pyfunction!(sign, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_cross, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_dist, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_len, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_norm, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_project, &m)?)?;
    m.add_function(wrap_pyfunction!(vec_psd, &m)?)?;
    march::register_module(&m)?;
    matrix::register_module(&m)?;
    optim::register_module(&m)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
