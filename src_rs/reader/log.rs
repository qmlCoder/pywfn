#![allow(non_snake_case)]

use std::collections::HashMap;

use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use crate::reader::{GauJob, Nics, Vibe};
use rswfn::reader::Reader;

#[pyclass]
pub struct LogReader {
    core: rswfn::reader::LogReader,
}

#[pymethods]
impl LogReader {
    /// 创建一个新的 LogReader 实例
    ///
    /// ## 参数
    /// - `path`: .log / .out 文件的路径
    ///
    /// ## 返回
    /// 返回一个 LogReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::LogReader::new(&path),
        };
        Ok(reader)
    }

    /// 读取作业关键词
    ///
    /// ## 返回
    /// 返回路由行内容
    pub fn read_job(&self) -> String {
        self.core.read_job().to_owned()
    }

    /// 解析作业关键词，得到任务类型、泛函与基组
    ///
    /// ## 返回
    /// 解析成功时返回 `GauJob` 对象，否则返回 None
    pub fn read_gau_job(&self) -> Option<GauJob> {
        self.core.read_gau_job().map(|core| GauJob { core })
    }

    /// 读取系数矩阵
    ///
    /// ## 参数
    /// - `spin`: 自旋类型，如 "Alpha" 或 "Beta"
    ///
    /// ## 返回
    /// 返回 (基函数所属原子索引, 基函数所属壳层索引, 基函数符号,
    /// 轨道能量列表, 轨道相关的字符串列表, 轨道系数矩阵)
    pub fn _read_coefs(
        &self,
        py: Python,
        spin: String,
    ) -> PyResult<(
        Vec<usize>,
        Vec<usize>,
        Vec<String>,
        Vec<f64>,
        Vec<String>,
        Py<PyArray2<f64>>,
    )> {
        match self.core._read_coefs(&spin) {
            Ok((ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat)) => Ok((
                ato_atms,
                ato_shls,
                ato_syms,
                obt_engs,
                obt_occs,
                cmat.into_pyarray(py).unbind(),
            )),
            Err(err) => Err(pyo3::exceptions::PyValueError::new_err(err)),
        }
    }

    /// 读取轨道系数（缓存结果）
    ///
    /// ## 返回
    /// 若文件中存在轨道系数则返回 (基函数所属原子索引, 基函数所属壳层索引,
    /// 基函数符号, 轨道能量列表, 轨道相关的字符串列表, 轨道系数矩阵)，否则返回 None
    pub fn read_coefs(
        &self,
        py: Python,
    ) -> Option<(
        Vec<usize>,
        Vec<usize>,
        Vec<String>,
        Vec<f64>,
        Vec<String>,
        Py<PyArray2<f64>>,
    )> {
        self.core.read_coefs().as_ref().map(|vals| {
            (
                vals.0.clone(),
                vals.1.clone(),
                vals.2.clone(),
                vals.3.clone(),
                vals.4.clone(),
                vals.5.clone().into_pyarray(py).unbind(),
            )
        })
    }

    /// 读取分子几何结构
    ///
    /// ## 参数
    /// - `irow`: 可选，指定起始行号，为 None 时自动查找
    ///
    /// ## 返回
    /// 返回 (原子序数列表, 原子坐标列表)
    #[pyo3(signature = (irow=None))]
    pub fn read_geome(&self, irow: Option<usize>) -> (Vec<usize>, Vec<[f64; 3]>) {
        self.core.read_geome(irow)
    }

    /// 读取所有优化步的几何结构
    ///
    /// ## 返回
    /// 返回 `Atoms` 对象列表，每个元素对应一步结构
    pub fn read_geomes(&self) -> Vec<Atoms> {
        self.core
            .read_geomes()
            .into_iter()
            .map(|core| Atoms { core })
            .collect()
    }

    /// 读取热力学能量数据
    ///
    /// ## 返回
    /// 返回能量列表，若没有能量信息则返回 None
    pub fn read_energys(&self) -> Option<Vec<f64>> {
        self.core.read_energys()
    }

    /// 读取振动频率
    ///
    /// ## 返回
    /// 返回频率列表，若没有频率信息则返回 None
    pub fn read_freqs(&self) -> Option<Vec<f64>> {
        self.core.read_freqs()
    }

    /// 读取对称矩阵
    ///
    /// ## 参数
    /// - `row`: 读取数据的起始行号
    /// - `loc`: 读取数据的起始列号
    ///
    /// ## 返回
    /// 返回 NumPy 二维数组形式的矩阵
    pub fn read_symMat(&self, py: Python, row: usize, loc: usize) -> Py<PyArray2<f64>> {
        self.core.read_symMat(row, loc).into_pyarray(py).unbind()
    }

    /// 读取密度矩阵
    ///
    /// ## 返回
    /// 返回 NumPy 二维数组形式的密度矩阵，若不存在则返回 None
    pub fn read_dmat(&self, py: Python) -> Option<Py<PyArray2<f64>>> {
        self.core
            .read_dmat()
            .map(|dmat| dmat.into_pyarray(py).unbind())
    }

    /// 读取重叠矩阵
    ///
    /// ## 返回
    /// 返回 NumPy 二维数组形式的重叠矩阵，若不存在则返回 None
    pub fn read_smat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_smat().map(|smat| smat.into_pyarray(py))
    }

    /// 读取动能矩阵
    ///
    /// ## 返回
    /// 返回 NumPy 二维数组形式的动能矩阵，若不存在则返回 None
    pub fn read_tmat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_tmat().map(|tmat| tmat.into_pyarray(py))
    }

    /// 读取势能矩阵
    ///
    /// ## 返回
    /// 返回 NumPy 二维数组形式的势能矩阵，若不存在则返回 None
    pub fn read_vmat<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.core.read_vmat().map(|vmat| vmat.into_pyarray(py))
    }

    /// 读取赝势基组的有效核电荷数
    ///
    /// ## 返回
    /// 返回 `{原子索引: 有效核电荷数}` 字典，若不存在则返回 None
    pub fn read_ecps(&self) -> Option<HashMap<usize, usize>> {
        self.core.read_ecps()
    }

    /// 读取振动模式信息
    ///
    /// ## 返回
    /// 返回 `Vibe` 对象列表，若文件不含振动信息则返回 None
    pub fn read_vibrations(&self) -> Option<Vec<Vibe>> {
        self.core
            .read_vibrations()
            .map(|vibes| vibes.into_iter().map(|core| Vibe { core }).collect())
    }

    /// 读取 NICS 数据
    ///
    /// ## 返回
    /// 返回 `Nics` 对象列表，若文件不含 NICS 信息则返回 None
    pub fn read_nics(&self) -> Option<Vec<Nics>> {
        self.core
            .read_nics()
            .map(|nics| nics.into_iter().map(|core| Nics { core }).collect())
    }

    /// 判断输出文件是否正常结束
    ///
    /// ## 返回
    /// 正常结束时返回 True，否则返回 False
    pub fn read_normalend(&self) -> bool {
        self.core.read_normalend()
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
}
