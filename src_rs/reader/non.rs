#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::Atoms;
use rswfn::reader::Reader;

#[pyclass]
pub struct NonReader {
    core: rswfn::reader::NonReader,
}

#[pymethods]
impl NonReader {
    /// 创建一个新的 NonReader 实例
    ///
    /// ## 参数
    /// - `path`: 文件路径
    ///
    /// ## 返回
    /// 返回一个 NonReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let result = Self {
            core: rswfn::reader::NonReader::new(&path),
        };
        Ok(result)
    }

    /// 设置原子索引和坐标
    ///
    /// ## 参数
    /// - `atms`: 原子索引列表
    /// - `xyzs`: 原子坐标列表，每个元素为 [x, y, z] 三元组
    pub fn set_atoms(&mut self, atms: Vec<usize>, xyzs: Vec<[f64; 3]>) {
        self.core.set_atoms(&atms, &xyzs);
    }
}
