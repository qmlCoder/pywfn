#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::{Atoms, Basis, Coefs};
use rswfn::reader::Reader;

/// Gaussian 任务信息（从 # 行解析得到）
#[pyclass]
pub struct GauJob {
    pub(crate) core: rswfn::reader::GauJob,
}

#[pymethods]
impl GauJob {
    /// 任务类型（如 opt、freq、sp 等）
    #[getter]
    pub fn job_type(&self) -> String {
        self.core.job_type.clone()
    }

    /// 泛函名称（如 B3LYP），未指定时为 None
    #[getter]
    pub fn functional(&self) -> Option<String> {
        self.core.functional.clone()
    }

    /// 基组名称（如 6-31G(d)），未指定时为 None
    #[getter]
    pub fn basis_set(&self) -> Option<String> {
        self.core.basis_set.clone()
    }
}

/// 解析高斯作业关键词行（形如 `#p b3lyp/6-31g(d) opt freq`）
///
/// ## 参数
/// - `line`: Gaussian 输入文件中的路由行（以 # 开头的那一行）
///
/// ## 返回
/// 解析成功时返回 `GauJob` 对象，否则返回 None
#[pyfunction]
pub fn parse_job(line: String) -> Option<GauJob> {
    rswfn::reader::parse_job(&line).map(|core| GauJob { core })
}

#[pyclass]
pub struct GjfReader {
    core: rswfn::reader::GjfReader,
}

#[pymethods]
impl GjfReader {
    /// 创建一个新的 GjfReader 实例
    ///
    /// ## 参数
    /// - `path`: .gjf 文件的路径
    ///
    /// ## 返回
    /// 返回一个 GjfReader 对象
    #[new]
    pub fn new(path: String) -> PyResult<Self> {
        let reader = Self {
            core: rswfn::reader::GjfReader::new(&path),
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

    /// 读取电荷和自旋多重度
    ///
    /// ## 返回
    /// 返回 (电荷, 自旋多重度)
    pub fn read_charge(&self) -> (isize, usize) {
        self.core.read_charge()
    }

    /// 读取作业关键词
    ///
    /// ## 返回
    /// 返回以 # 开头的路由行内容
    pub fn read_job(&self) -> String {
        self.core.read_job()
    }

    /// 解析作业关键词，得到任务类型、泛函与基组
    ///
    /// ## 返回
    /// 解析成功时返回 `GauJob` 对象，否则返回 None
    pub fn read_gau_job(&self) -> Option<GauJob> {
        self.core.read_gau_job().map(|core| GauJob { core })
    }
}
