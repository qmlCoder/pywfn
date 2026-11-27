#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use crate::reader::BasData;
use rswfn::reader::Reader;

#[pyclass]
pub struct MdeReader {
    core: rswfn::reader::MdeReader,
}

#[pymethods]
impl MdeReader {
    /// 创建一个新的 MdeReader 实例
    ///
    /// ## 参数
    /// - `path`: .mde 文件的路径
    ///
    /// ## 返回
    /// 返回一个 MdeReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::MdeReader::new(&path),
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
    /// 返回 (原子序数列表, 原子坐标列表)
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        self.core.read_geome()
    }

    /// 读取基组信息（带缓存）
    ///
    /// ## 返回
    /// 返回 (基函数所属原子索引, 基函数所属壳层索引, 基函数符号, 基函数原始数据列表)
    pub fn read_basis(&self) -> (Vec<usize>, Vec<usize>, Vec<String>, Vec<BasData>) {
        let (ato_atms, ato_shls, ato_syms, ato_bass) = self.core.read_basis();
        let ato_bass = ato_bass
            .iter()
            .map(|bas| BasData::from_core(*bas))
            .collect();
        (
            ato_atms.to_vec(),
            ato_shls.to_vec(),
            ato_syms.to_vec(),
            ato_bass,
        )
    }

    /// 读取轨道系数
    ///
    /// ## 返回
    /// 返回 (轨道能量列表, 轨道占据数列表, 轨道系数矩阵)
    pub fn read_coefs(&self) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
        self.core.read_coefs()
    }

    /// 读取电子数
    ///
    /// ## 返回
    /// 返回 `[α 电子数, β 电子数]`
    pub fn read_neles(&self) -> [usize; 2] {
        self.core.read_neles()
    }
}
