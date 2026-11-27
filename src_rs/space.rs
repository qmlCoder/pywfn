#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "space")?;
    m.add_function(wrap_pyfunction!(a2m_weits, &m)?)?;
    m.add_function(wrap_pyfunction!(calc_EIEBA, &m)?)?;
    m.add_function(wrap_pyfunction!(get_grids, &m)?)?;
    parent_module.add_submodule(&m)
}

/// 将原子分区的格点权重重新分配到整个分子
///
/// ## 参数
/// - `iatm`: 目标原子索引
/// - `atm_grids`: 原子分区的格点坐标
/// - `atm_weits`: 原子分区的格点权重
/// - `atm_pos`: 所有原子的位置
/// - `atm_rad`: 所有原子的半径
///
/// ## 返回
/// - 重新分配后的权重列表
#[pyfunction]
pub fn a2m_weits(
    iatm: usize,
    atm_grids: Vec<[f64; 3]>,
    atm_weits: Vec<f64>,
    atm_pos: Vec<[f64; 3]>,
    atm_rad: Vec<f64>,
) -> Vec<f64> {
    rswfn::space::a2m_weits(iatm, &atm_grids, &atm_weits, &atm_pos, &atm_rad)
}

/// 计算两个片段之间的 EIEBA 值
///
/// ## 参数
/// - `fragA`: 片段 A 的格点数据，每个元素为 [x, y, z, 值]
/// - `fragB`: 片段 B 的格点数据
/// - `rhoA`: 片段 A 的密度值
/// - `rhoB`: 片段 B 的密度值
/// - `gridsA`: 片段 A 的格点坐标
/// - `gridsB`: 片段 B 的格点坐标
/// - `weitsA`: 片段 A 的格点权重
/// - `weitsB`: 片段 B 的格点权重
///
/// ## 返回
/// - 三个计算结果组成的元组
#[pyfunction]
pub fn calc_EIEBA(
    fragA: Vec<[f64; 4]>,
    fragB: Vec<[f64; 4]>,
    rhoA: Vec<f64>,
    rhoB: Vec<f64>,
    gridsA: Vec<[f64; 3]>,
    gridsB: Vec<[f64; 3]>,
    weitsA: Vec<f64>,
    weitsB: Vec<f64>,
) -> (f64, f64, f64) {
    rswfn::space::calc_EIEBA(
        &fragA, &fragB, &rhoA, &rhoB, &gridsA, &gridsB, &weitsA, &weitsB,
    )
}

/// 生成三维格点坐标
///
/// ## 参数
/// - `nx`: x 方向格点数
/// - `ny`: y 方向格点数
/// - `nz`: z 方向格点数
///
/// ## 返回
/// - 格点坐标列表
#[pyfunction]
pub fn get_grids(nx: usize, ny: usize, nz: usize) -> Vec<[f64; 3]> {
    rswfn::space::get_grids(nx, ny, nz)
}
