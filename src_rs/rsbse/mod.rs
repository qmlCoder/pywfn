//! Basis Set Exchange（BSE）基组读取模块。
//!
//! 对应 `rswfn::rsbse`，负责从网络、文件或文本中读取基组并转换为 `Basis`。

use crate::base::{Atoms, Basis};
use pyo3::prelude::*;
use rswfn;

/// 从基组文件读取基组
///
/// 把 `path` 指定的基组文件内容解析为 `Basis`，仅保留 `atoms` 中出现的元素。
/// 该函数不设置超时，读取大文件或网络文件可能耗时较长。
///
/// ## 参数
/// - `path`: 基组文件路径（如 BSE 下载的 .gbs 文件）
/// - `atoms`: 原子集合对象，用于确定需要哪些元素的基组
///
/// ## 返回
/// 成功时返回 Basis 对象，解析失败时返回 None
#[pyfunction]
pub fn from_file(path: String, atoms: PyRef<'_, Atoms>) -> Option<Basis> {
    rswfn::rsbse::from_file(&path, &atoms.core).map(|basis| Basis { core: basis })
}

/// 按名称从 Basis Set Exchange 下载并解析基组
///
/// 名称通过 `rswfn::rsbse::BASIS_SETS` 中的别名表匹配。
/// 该函数会访问网络，且不设置超时，可能很慢甚至失败。
///
/// ## 参数
/// - `name`: 基组名称，如 "def2-SVP"
/// - `atoms`: 原子集合对象，用于确定需要哪些元素的基组
///
/// ## 返回
/// 成功时返回 Basis 对象，名称不存在或下载失败时返回 None
#[pyfunction]
pub fn from_name(name: String, atoms: PyRef<'_, Atoms>) -> Option<Basis> {
    rswfn::rsbse::from_name(&name, &atoms.core).map(|basis| Basis { core: basis })
}

/// 从基组文本内容解析基组
///
/// ## 参数
/// - `text`: 基组文本内容（Gaussian 格式）
/// - `atoms`: 原子集合对象，用于确定需要哪些元素的基组
///
/// ## 返回
/// 成功时返回 Basis 对象，解析失败时返回 None
#[pyfunction]
pub fn from_text(text: String, atoms: PyRef<'_, Atoms>) -> Option<Basis> {
    rswfn::rsbse::from_text(&text, &atoms.core).map(|basis| Basis { core: basis })
}

/// 从基组文本的逐行内容中解析壳层数据
///
/// ## 参数
/// - `lines`: 基组文本按行拆分后的列表
///
/// ## 返回
/// 成功时返回列表，每个元素为 (原子序号, 壳层序号, 角动量, 高斯指数, 收缩系数)；
/// 解析失败时返回 None
#[pyfunction]
pub fn read_basis(lines: Vec<String>) -> Option<Vec<(usize, usize, usize, f64, f64)>> {
    rswfn::rsbse::read_basis(&lines)
}

/// 返回所有支持的基组名称
///
/// ## 返回
/// 返回基组名称列表，可直接传给 `from_name`
#[pyfunction]
pub fn supported_names() -> Vec<String> {
    rswfn::rsbse::supported_names()
}

/// 返回内置的基组别名表
///
/// ## 返回
/// 返回 (别名, BSE 名称) 二元组列表
fn basis_sets() -> Vec<(&'static str, &'static str)> {
    rswfn::rsbse::BASIS_SETS.to_vec()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsbse")?;
    m.add_function(wrap_pyfunction!(from_file, &m)?)?;
    m.add_function(wrap_pyfunction!(from_name, &m)?)?;
    m.add_function(wrap_pyfunction!(from_text, &m)?)?;
    m.add_function(wrap_pyfunction!(read_basis, &m)?)?;
    m.add_function(wrap_pyfunction!(supported_names, &m)?)?;
    // 静态表以模块属性的形式暴露，对应 rswfn::rsbse::BASIS_SETS
    m.add("BASIS_SETS", basis_sets())?;
    parent_module.add_submodule(&m)
}
