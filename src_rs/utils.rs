#![allow(non_snake_case)]

use ndarray::Array2;
use numpy::{PyArray2, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "utils")?;
    m.add_function(wrap_pyfunction!(rswfn_version, &m)?)?;
    m.add_function(wrap_pyfunction!(all_same, &m)?)?;
    m.add_function(wrap_pyfunction!(comp_matrix, &m)?)?;
    m.add_function(wrap_pyfunction!(file_name, &m)?)?;
    m.add_function(wrap_pyfunction!(file_type, &m)?)?;
    m.add_function(wrap_pyfunction!(fort_form, &m)?)?;
    m.add_function(wrap_pyfunction!(has_nan, &m)?)?;
    m.add_function(wrap_pyfunction!(posi_matrix, &m)?)?;
    m.add_function(wrap_pyfunction!(same_matrix, &m)?)?;
    m.add_function(wrap_pyfunction!(show_matrix, &m)?)?;
    m.add_function(wrap_pyfunction!(symm_matrix, &m)?)?;
    parent_module.add_submodule(&m)
}

/// 返回底层依赖的 rswfn 库版本号。
///
/// 该值在构建时从 Cargo.lock 中读取，反映实际编译链接的 rswfn 版本。
#[pyfunction]
pub fn rswfn_version() -> &'static str {
    env!("RSWFN_VERSION")
}

/// 逐个比较两个列表对应位置的值是否完全相同。
///
/// ## 参数
/// - `vals0`: 第一个列表
/// - `vals1`: 第二个列表
///
/// ## 返回
/// - 长度相同且每个元素都相等时返回 `True`，否则返回 `False`
#[pyfunction]
pub fn all_same(vals0: Vec<f64>, vals1: Vec<f64>) -> bool {
    rswfn::utils::all_same(&vals0, &vals1)
}

/// 比较两个矩阵，打印首个不相同的元素位置。
///
/// ## 参数
/// - `mat1`: 第一个矩阵
/// - `mat2`: 第二个矩阵
/// - `show`: 最大打印个数，非正数表示不打印
#[pyfunction]
pub fn comp_matrix(mat1: Bound<'_, PyArray2<f64>>, mat2: Bound<'_, PyArray2<f64>>, show: isize) {
    let mat1 = mat1.readonly().as_array().to_owned();
    let mat2 = mat2.readonly().as_array().to_owned();
    rswfn::utils::comp_matrix(&mat1, &mat2, show)
}

/// 获取文件路径中的文件名（不含目录部分）。
///
/// ## 参数
/// - `path`: 文件路径
#[pyfunction]
pub fn file_name(path: &str) -> String {
    rswfn::utils::file_name(path)
}

/// 获取文件路径的扩展名。
///
/// ## 参数
/// - `path`: 文件路径
#[pyfunction]
pub fn file_type(path: &str) -> String {
    rswfn::utils::file_type(path)
}

/// 按 Fortran 风格格式化浮点数。
///
/// ## 参数
/// - `val`: 需要格式化的值
/// - `pec`: 小数位数
#[pyfunction]
pub fn fort_form(val: f64, pec: usize) -> String {
    rswfn::utils::fort_form(val, pec)
}

/// 判断矩阵中是否存在 `NaN`。
///
/// ## 参数
/// - `matrix`: 待检查的矩阵
#[pyfunction]
pub fn has_nan(matrix: Bound<'_, PyArray2<f64>>) -> bool {
    let matrix = matrix.readonly().as_array().to_owned();
    rswfn::utils::has_nan(&matrix)
}

/// 判断矩阵是否所有元素都为正。
///
/// ## 参数
/// - `mat`: 待检查的矩阵
#[pyfunction]
pub fn posi_matrix(mat: Bound<'_, PyArray2<f64>>) -> bool {
    let mat = mat.readonly().as_array().to_owned();
    rswfn::utils::posi_matrix(&mat)
}

/// 判断两个矩阵是否完全相同。
///
/// ## 参数
/// - `mat1`: 第一个矩阵
/// - `mat2`: 第二个矩阵
#[pyfunction]
pub fn same_matrix(mat1: Bound<'_, PyArray2<f64>>, mat2: Bound<'_, PyArray2<f64>>) -> bool {
    let mat1 = mat1.readonly().as_array().to_owned();
    let mat2 = mat2.readonly().as_array().to_owned();
    rswfn::utils::same_matrix(&mat1, &mat2)
}

/// 打印矩阵内容。
///
/// ## 参数
/// - `mat`: 需要打印的矩阵
#[pyfunction]
pub fn show_matrix(mat: Bound<'_, PyArray2<f64>>) {
    let mat = mat.readonly().as_array().to_owned();
    rswfn::utils::show_matrix(&mat)
}

/// 判断矩阵是否对称。
///
/// ## 参数
/// - `mat`: 待检查的矩阵
#[pyfunction]
pub fn symm_matrix(mat: Bound<'_, PyArray2<f64>>) -> bool {
    let mat = mat.readonly().as_array().to_owned();
    rswfn::utils::symm_matrix(&mat)
}

pub fn arr2vec(arr: &Array2<f64>) -> Vec<Vec<f64>> {
    let nrow = arr.nrows();
    let ncol = arr.ncols();
    let mut vec = vec![vec![0.0; ncol]; nrow];
    for i in 0..nrow {
        for j in 0..ncol {
            vec[i][j] = arr[[i, j]];
        }
    }
    vec
}
