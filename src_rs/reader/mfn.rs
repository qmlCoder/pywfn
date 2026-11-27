#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use crate::reader::BasData;
use rswfn::reader::Reader;

/// Multiwfn 波函数文件(.mfn)读取器
#[pyclass]
pub struct MfnReader {
    core: rswfn::reader::MfnReader,
}

#[pymethods]
impl MfnReader {
    /// 创建一个新的 MfnReader 实例
    ///
    /// ## 参数
    /// - `path`: .mfn 文件的路径
    ///
    /// ## 返回
    /// 返回一个 MfnReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::MfnReader::new(&path),
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
    /// 返回原子序数、笛卡尔坐标（内部使用 Bohr）和有效核电荷数。
    ///
    /// ## 返回
    /// 返回 (原子序数列表, 原子坐标列表, 有效核电荷数列表)
    pub fn read_geome(&self) -> (Vec<usize>, Vec<[f64; 3]>, Vec<f64>) {
        self.core.read_geome()
    }

    /// 读取基函数信息
    ///
    /// 返回每个原子轨道对应的原子索引、壳层索引、符号，原始基组数据
    /// 以及每个壳层的基函数形式（"car" 或 "sph"）。
    ///
    /// ## 返回
    /// 返回 (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
    /// 基函数原始数据列表, 基函数形式列表)
    pub fn read_basis(
        &self,
    ) -> (
        Vec<usize>,
        Vec<usize>,
        Vec<String>,
        Vec<BasData>,
        Vec<String>,
    ) {
        let (ato_atms, ato_shls, ato_syms, ato_bass, ato_forms) = self.core.read_basis();
        let ato_bass = ato_bass
            .into_iter()
            .map(|core| BasData::from_core(core))
            .collect();
        (ato_atms, ato_shls, ato_syms, ato_bass, ato_forms)
    }

    /// 读取轨道系数、能量和占据数
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

    /// 读取 α 和 β 电子数
    ///
    /// ## 返回
    /// 返回 `[α 电子数, β 电子数]`
    pub fn read_neles(&self) -> [usize; 2] {
        self.core.read_neles()
    }
}
