#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use rswfn::reader::Reader;

#[pyclass]
pub struct CubReader {
    core: rswfn::reader::CubReader,
}

#[pymethods]
impl CubReader {
    /// 创建一个新的 CubReader 实例
    ///
    /// ## 参数
    /// - `path`: .cub 文件的路径
    ///
    /// ## 返回
    /// 返回一个 CubReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::CubReader::new(&path),
        };
        Ok(reader)
    }

    /// 读取分子中的原子信息
    ///
    /// ## 返回
    /// 返回 `Atoms` 对象，包含原子符号、坐标等信息
    pub fn get_atoms(&self) -> Atoms {
        let atoms = self.core.get_atoms();
        Atoms { core: atoms }
    }

    /// 读取基组信息（如存在）
    ///
    /// ## 返回
    /// 返回 `Option<Basis>`，若文件中包含基组信息则返回 Some(Basis)，否则返回 None
    pub fn get_basis(&self) -> Option<Basis> {
        let basis = self.core.get_basis();
        if let Some(basis) = basis {
            Some(Basis { core: basis })
        } else {
            None
        }
    }

    /// 读取分子轨道系数（如存在）
    ///
    /// ## 返回
    /// 返回 `Option<Coefs>`，若文件中包含轨道系数则返回 Some(Coefs)，否则返回 None
    pub fn get_coefs(&self) -> Option<Coefs> {
        let coefs = self.core.get_coefs();
        if let Some(coefs) = coefs {
            Some(Coefs { core: coefs })
        } else {
            None
        }
    }

    /// 读取格点信息
    ///
    /// ## 返回
    /// 返回 `CubeGrid` 对象，包含三个方向的格点数、起点坐标和步长
    pub fn read_cube(&self) -> crate::props::gridprop::CubeGrid {
        crate::props::gridprop::CubeGrid::from_core(self.core.read_cube())
    }

    /// 读取原子信息
    ///
    /// ## 返回
    /// 返回 (原子序数列表, 原子电荷列表, 原子坐标列表)
    pub fn read_atms(&self) -> (Vec<usize>, Vec<f64>, Vec<[f64; 3]>) {
        self.core.read_atms()
    }

    /// 读取格点数值
    ///
    /// ## 返回
    /// 返回 (轨道索引列表, 每个轨道在全部格点上的数值列表)
    pub fn read_vals(&self) -> (Vec<usize>, Vec<Vec<f64>>) {
        self.core.read_vals()
    }
}
