#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

/// 使用 Marching Cubes 算法提取等值面
///
/// 从三维标量场中提取给定等值面的三角形网格。
///
/// ## 参数
/// - `shape`: 三个方向上的格点数量 [nx, ny, nz]
/// - `grids`: 所有格点的坐标列表
/// - `values`: 所有格点上的标量值列表
/// - `isov`: 等值面的值
///
/// ## 返回
/// 返回一个元组 (vertices, triangles)，其中
/// - `vertices`: 三角形顶点坐标列表
/// - `triangles`: 三角形顶点索引列表（每三个为一组）
#[pyfunction]
pub fn march_cube_v1(
    shape: [usize; 3],
    grids: Vec<[f64; 3]>,
    values: Vec<f64>,
    isov: f64,
) -> PyResult<(Vec<[f64; 3]>, Vec<isize>)> {
    let res = rswfn::march::mc::march_cube_v1(&shape, &grids, &values, isov);
    Ok(res)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "march")?; // 为当前的rs文件创建一个子模块
    m.add_function(wrap_pyfunction!(march_cube_v1, &m)?)?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
