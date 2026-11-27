#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use rswfn::reader::Reader;

#[pyclass]
pub struct DmoReader {
    core: rswfn::reader::DmoReader,
}

#[pymethods]
impl DmoReader {
    /// 创建一个新的 DmoReader 实例
    ///
    /// ## 参数
    /// - `path`: .dmo 文件的路径
    ///
    /// ## 返回
    /// 返回一个 DmoReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::DmoReader::new(&path),
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

    /// 读取原子信息（带缓存）
    ///
    /// ## 返回
    /// 返回 (原子序数列表, 原子坐标列表)
    pub fn read_atms(&self) -> (Vec<usize>, Vec<[f64; 3]>) {
        (*self.core.read_atms()).clone()
    }

    /// 读取原子轨道信息（带缓存）
    ///
    /// ## 返回
    /// 返回每个原子轨道的 (原子索引, 壳层索引, 符号, 坐标)
    pub fn read_atos(&self) -> Vec<(usize, usize, String, [f64; 3])> {
        (*self.core.read_atos()).clone()
    }

    /// 读取键连矩阵
    ///
    /// ## 返回
    /// 返回键连矩阵（NumPy 二维数组）
    pub fn read_BM(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.read_BM().into_pyarray(py).unbind()
    }

    /// 读取轨道系数和能量
    ///
    /// ## 返回
    /// 返回 (轨道能量列表, 轨道系数矩阵)
    pub fn read_coefs(&self, py: Python) -> (Vec<f64>, Py<PyArray2<f64>>) {
        let (obt_engs, cmat) = self.core.read_coefs();
        (obt_engs, cmat.into_pyarray(py).unbind())
    }
}
