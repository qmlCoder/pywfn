#![allow(non_snake_case)]

use std::collections::HashMap;

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use rswfn::reader::Reader;

#[pyclass]
pub struct OutReader {
    core: rswfn::reader::OutReader,
}

#[pymethods]
impl OutReader {
    /// 创建一个新的 OutReader 实例
    ///
    /// ## 参数
    /// - `path`: .out 文件的路径
    ///
    /// ## 返回
    /// 返回一个 OutReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::OutReader::new(&path),
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

    /// 读取基组信息
    ///
    /// ## 返回
    /// 返回 `Basis` 对象，包含基函数、收缩系数等信息
    pub fn get_basis(&self) -> Option<Basis> {
        let basis = self.core.get_basis();
        if let Some(basis) = basis {
            Some(Basis { core: basis })
        } else {
            None
        }
    }

    /// 读取分子轨道系数
    ///
    /// ## 返回
    /// 返回 `Coefs` 对象，包含分子轨道系数矩阵
    pub fn get_coefs(&self) -> Option<Coefs> {
        let coefs = self.core.get_coefs();
        if let Some(coefs) = coefs {
            Some(Coefs { core: coefs })
        } else {
            None
        }
    }

    /// 读取几何结构
    ///
    /// ## 返回
    /// 返回 (原子符号列表, 原子坐标列表)
    pub fn read_geome(&self) -> (Vec<String>, Vec<[f64; 3]>) {
        self.core.read_geome()
    }

    /// 读取基组信息
    ///
    /// ## 返回
    /// 返回 `{原子索引: [(角动量, 壳层索引, 指数, 系数), ...]}` 字典
    pub fn read_basis(&self) -> HashMap<usize, Vec<(i32, usize, f64, f64)>> {
        self.core.read_basis()
    }

    /// 读取分子轨道系数
    ///
    /// ## 返回
    /// 返回 (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
    /// 轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
    pub fn read_coefs(
        &self,
        py: Python,
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<String>,
        Vec<f64>,
        Vec<f64>,
        Py<PyArray2<f64>>,
    ) {
        let (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) = self.core.read_coefs();
        (
            ato_atms,
            ato_shls,
            ato_syms,
            obt_engs,
            obt_occs,
            cmat.into_pyarray(py).unbind(),
        )
    }

    /// 读取轨道能量和占据数
    ///
    /// ## 返回
    /// 返回 (轨道能量列表, 轨道占据数列表)
    pub fn read_obtEngs(&self) -> (Vec<f64>, Vec<f64>) {
        self.core.read_obtEngs()
    }
}
