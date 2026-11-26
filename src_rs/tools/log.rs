use pyo3::prelude::*;

use rswfn;

#[pyclass]
struct Tool {
    core: rswfn::tools::log::Tool,
}

#[pymethods]
impl Tool {
    /// 创建一个新的 Tool 实例
    ///
    /// ## 返回
    /// 返回一个 Tool 对象
    #[new]
    fn new() -> Tool {
        Tool {
            core: rswfn::tools::log::Tool::new(),
        }
    }

    /// 从高斯 .log 文件中读取 SI（稳定性指标）信息
    ///
    /// ## 参数
    /// - `paths`: .log 文件路径列表
    ///
    /// ## 返回
    /// 返回 SI 字符串，若未找到则返回空字符串
    pub fn get_SI(&self, paths: Vec<String>) -> String {
        self.core.get_SI(&paths)
    }

    /// 拆分优化轨迹，提取每一步的原子索引、坐标和能量
    ///
    /// ## 参数
    /// - `path`: .log 文件路径
    ///
    /// ## 返回
    /// 返回优化过程中每一步的 (原子索引, 坐标, 能量) 元组列表
    pub fn split_opt(&self, path: String) -> Vec<(Vec<usize>, Vec<[f64; 3]>, f64)> {
        self.core.split_opt(path)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "log")?; // 为当前的rs文件创建一个子模块
    m.add_class::<Tool>()?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
