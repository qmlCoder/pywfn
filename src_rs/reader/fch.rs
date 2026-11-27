#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs, Stm};
use rswfn::reader::Reader;

#[pyclass]
pub struct FchReader {
    core: rswfn::reader::FchReader,
}

#[pymethods]
impl FchReader {
    /// 创建一个新的 FchReader 实例
    ///
    /// ## 参数
    /// - `path`: .fch 文件的路径
    ///
    /// ## 返回
    /// 返回一个 FchReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::FchReader::new(&path),
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

    /// 读取原子序数
    ///
    /// ## 返回
    /// 返回每个原子的原子序数列表（带符号，用于标记赝势等特殊原子）
    pub fn read_atms(&self) -> Vec<isize> {
        self.core.read_atms()
    }

    /// 读取原子核电荷（无 ECP 时为原子序数）
    ///
    /// ## 返回
    /// 返回每个原子的电荷列表
    pub fn read_crgs(&self) -> Vec<f64> {
        self.core.read_crgs()
    }

    /// 读取笛卡尔坐标
    ///
    /// ## 返回
    /// 返回每个原子的 [x, y, z] 坐标列表
    pub fn read_xyzs(&self) -> Vec<[f64; 3]> {
        self.core.read_xyzs()
    }

    /// 读取局部坐标系
    ///
    /// ## 返回
    /// 若文件中存在局部坐标系信息则返回 `Stm` 列表，否则返回 None
    pub fn read_stms(&self) -> Option<Vec<Stm>> {
        self.core
            .read_stms()
            .map(|stms| stms.into_iter().map(|core| Stm { core }).collect())
    }

    /// 读取轨道系数及其他轨道相关信息
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
}
