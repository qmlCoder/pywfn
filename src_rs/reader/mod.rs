mod cub;
mod dmo;
mod fch;
mod gjf;
mod inp;
mod log;
mod mde;
mod mfn;
mod mol;
mod non;
mod out;
mod sdf;
mod xyz;

use std::collections::HashMap;

use pyo3::prelude::*;

pub use cub::CubReader;
pub use dmo::DmoReader;
pub use fch::FchReader;
pub use gjf::{parse_job, GauJob, GjfReader};
pub use inp::InpReader;
pub use log::LogReader;
pub use mde::MdeReader;
pub use mfn::MfnReader;
pub use mol::MolReader;
pub use non::NonReader;
pub use out::OutReader;
pub use sdf::SdfReader;
pub use xyz::XyzReader;

use crate::base::Basis;
use rswfn;

/// 支持读取的文件类型
///
/// 与 rswfn 的 `reader::FileType` 枚举一一对应
#[pyclass(eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FileType {
    /// Gaussian 输入文件 (.gjf)
    Gjf,
    /// Gaussian 输出文件 (.log / .out)
    Log,
    /// Gaussian 格式化检查点文件 (.fch / .fchk)
    Fch,
    /// MDL 分子文件 (.mol)
    Mol,
    /// MDL 结构数据文件 (.sdf)
    Sdf,
    /// 标准 XYZ 坐标文件 (.xyz)
    Xyz,
    /// Molden 波函数文件 (.molden)
    Mde,
    /// Multiwfn 波函数文件 (.mfn)
    Mfn,
    /// DMol3 输出文件 (.dmo / .dmi)
    Dmo,
    /// 通用输入文件 (.inp)
    Inp,
    /// 立方体文件 (.cub)
    Cub,
}

/// 把 rswfn 的文件类型枚举转换为 pywfn 的枚举
///
/// ## 参数
/// - `filetype`: rswfn 中的文件类型
///
/// ## 返回
/// 返回 pywfn 中对应的 `FileType`
#[allow(unreachable_patterns)]
fn filetype_from_core(filetype: rswfn::reader::FileType) -> FileType {
    match filetype {
        rswfn::reader::FileType::Gjf => FileType::Gjf,
        rswfn::reader::FileType::Log => FileType::Log,
        rswfn::reader::FileType::Fch => FileType::Fch,
        rswfn::reader::FileType::Mol => FileType::Mol,
        rswfn::reader::FileType::Sdf => FileType::Sdf,
        rswfn::reader::FileType::Xyz => FileType::Xyz,
        rswfn::reader::FileType::Mde => FileType::Mde,
        rswfn::reader::FileType::Mfn => FileType::Mfn,
        rswfn::reader::FileType::Dmo => FileType::Dmo,
        rswfn::reader::FileType::Inp => FileType::Inp,
        rswfn::reader::FileType::Cub => FileType::Cub,
        // 上游若新增变体则兜底为 Log，避免此处编译失败
        _ => FileType::Log,
    }
}

/// 基函数的原始数据
///
/// 对应 rswfn 中的类型别名 `BasData = (usize, usize, usize, f64, f64)`，
/// 即 (原子索引, 壳层索引, 角动量, 高斯指数, 收缩系数)
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct BasData {
    /// 所属原子索引
    pub atm: usize,
    /// 所属壳层索引
    pub shl: usize,
    /// 角动量量子数
    pub ang: usize,
    /// 高斯指数
    pub alp: f64,
    /// 收缩系数
    pub coe: f64,
}

impl BasData {
    /// 由 rswfn 的元组数据构造
    ///
    /// ## 参数
    /// - `core`: rswfn 中的 `(atm, shl, ang, alp, coe)` 元组
    ///
    /// ## 返回
    /// 返回 `BasData` 对象
    pub(crate) fn from_core(core: rswfn::reader::BasData) -> Self {
        Self {
            atm: core.0,
            shl: core.1,
            ang: core.2,
            alp: core.3,
            coe: core.4,
        }
    }

    /// 转换为 rswfn 的元组数据
    ///
    /// ## 返回
    /// 返回 `(atm, shl, ang, alp, coe)` 元组
    pub(crate) fn to_core(&self) -> rswfn::reader::BasData {
        (self.atm, self.shl, self.ang, self.alp, self.coe)
    }
}

#[pymethods]
impl BasData {
    /// 创建一个基函数原始数据对象
    ///
    /// ## 参数
    /// - `atm`: 所属原子索引
    /// - `shl`: 所属壳层索引
    /// - `ang`: 角动量量子数
    /// - `alp`: 高斯指数
    /// - `coe`: 收缩系数
    #[new]
    pub fn new(atm: usize, shl: usize, ang: usize, alp: f64, coe: f64) -> Self {
        Self {
            atm,
            shl,
            ang,
            alp,
            coe,
        }
    }

    /// 所属原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.atm
    }

    /// 所属壳层索引
    #[getter]
    pub fn shl(&self) -> usize {
        self.shl
    }

    /// 角动量量子数
    #[getter]
    pub fn ang(&self) -> usize {
        self.ang
    }

    /// 高斯指数
    #[getter]
    pub fn alp(&self) -> f64 {
        self.alp
    }

    /// 收缩系数
    #[getter]
    pub fn coe(&self) -> f64 {
        self.coe
    }

    /// 转换为 Python 元组
    ///
    /// ## 返回
    /// 返回 `(atm, shl, ang, alp, coe)` 元组
    pub fn as_tuple(&self) -> (usize, usize, usize, f64, f64) {
        (self.atm, self.shl, self.ang, self.alp, self.coe)
    }

    pub fn __repr__(&self) -> String {
        format!(
            "BasData(atm={}, shl={}, ang={}, alp={}, coe={})",
            self.atm, self.shl, self.ang, self.alp, self.coe
        )
    }
}

/// 文件中的关键字标记表
///
/// 记录每个关键字（如 "NAtoms"）在文件行列表中的行号，
/// 用于快速定位数据块。
#[pyclass]
pub struct Marks {
    core: rswfn::reader::Marks,
}

#[pymethods]
impl Marks {
    /// 创建新的 Marks 实例
    ///
    /// ## 返回
    /// 返回一个空的 `Marks` 对象
    #[new]
    pub fn new() -> Self {
        Self {
            core: rswfn::reader::Marks::new(),
        }
    }

    /// 获取全部标记及其行号
    ///
    /// ## 返回
    /// 返回 `{标记名: 行号}` 字典
    #[getter]
    pub fn marks(&self) -> HashMap<String, usize> {
        self.core.marks.clone()
    }

    /// 添加一个标记
    ///
    /// ## 参数
    /// - `marker`: 标记名称
    /// - `number`: 标记所在的行号
    pub fn add(&mut self, marker: String, number: usize) {
        self.core.add(marker, number);
    }
}

#[pyclass]
pub struct Nics {
    core: rswfn::reader::Nics,
}

#[pymethods]
impl Nics {
    #[getter]
    fn iso(&self) -> f64 {
        self.core.iso
    }

    #[getter]
    fn ani(&self) -> f64 {
        self.core.ani
    }

    #[getter]
    fn mat(&self) -> [f64; 9] {
        self.core.mat
    }

    #[getter]
    fn eig(&self) -> [f64; 3] {
        self.core.eig
    }
}

/// 一个正则振动模式
///
/// 包含振动频率、约化质量、力常数、红外强度以及各原子的振动位移
#[pyclass]
pub struct Vibe {
    core: rswfn::reader::Vibe,
}

#[pymethods]
impl Vibe {
    /// 振动频率
    #[getter]
    pub fn freq(&self) -> f64 {
        self.core.freq
    }

    /// 约化质量
    #[getter]
    pub fn redm(&self) -> f64 {
        self.core.redm
    }

    /// 力常数
    #[getter]
    pub fn frcc(&self) -> f64 {
        self.core.frcc
    }

    /// 红外强度
    #[getter]
    pub fn irin(&self) -> f64 {
        self.core.irin
    }

    /// 每个原子的振动位移向量
    #[getter]
    pub fn vibs(&self) -> Vec<[f64; 3]> {
        self.core.vibs.clone()
    }

    /// 获取振动幅度最大的原子及其位移方向
    ///
    /// ## 返回
    /// 返回 (原子索引, 该原子的位移向量, 该原子的位移大小)
    pub fn get_max(&self) -> (usize, [f64; 3], f64) {
        self.core.get_max()
    }
}

/// 读取文件的全部行
///
/// 打开指定路径的文件，逐行读取并返回字符串列表。
///
/// ## 参数
/// - `path`: 文件路径
///
/// ## 返回
/// 返回文件每一行组成的字符串列表
#[pyfunction]
pub fn read_lines(path: String) -> Vec<String> {
    rswfn::reader::read_lines(&path)
}

/// 根据电荷和自旋多重度计算 α 和 β 电子数
///
/// ## 参数
/// - `charge`: 体系电荷
/// - `multip`: 自旋多重度
/// - `total`: 总电子数
///
/// ## 返回
/// 返回 `[α 电子数, β 电子数]`
#[pyfunction]
pub fn charge_number(charge: isize, multip: usize, total: usize) -> [usize; 2] {
    rswfn::reader::charge_number(charge, multip, total)
}

/// 判断文件类型
///
/// 根据文件扩展名返回对应的 `FileType` 枚举值。
///
/// ## 参数
/// - `path`: 文件路径
///
/// ## 返回
/// 返回对应的 `FileType` 枚举值
#[pyfunction]
pub fn get_fileType(path: String) -> FileType {
    filetype_from_core(rswfn::reader::get_fileType(&path))
}

/// 从原始基组数据构建 `Basis` 对象
///
/// ## 参数
/// - `datas`: 每个基函数的原始数据（原子索引、壳层索引、角量子数、指数、收缩系数）
/// - `forms`: 每个壳层的基函数形式列表
///
/// ## 返回
/// 返回 `Basis` 基组对象
#[pyfunction]
pub fn shl_basis(datas: Vec<BasData>, forms: Vec<String>) -> Basis {
    let cores: Vec<rswfn::reader::BasData> = datas.iter().map(|data| data.to_core()).collect();
    let basis = rswfn::reader::shl_basis(&cores, &forms);
    Basis { core: basis }
}

/// 确定每个壳层的基函数形式（球谐或笛卡尔）
///
/// ## 参数
/// - `ato_atms`: 每个基函数所属的原子索引
/// - `ato_shls`: 每个基函数所属的壳层索引
/// - `ato_syms`: 每个基函数的符号
///
/// ## 返回
/// 返回每个基函数对应的形式字符串列表
#[pyfunction]
pub fn shl_forms(ato_atms: Vec<usize>, ato_shls: Vec<usize>, ato_syms: Vec<String>) -> Vec<String> {
    rswfn::reader::shl_forms(&ato_atms, &ato_shls, &ato_syms)
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "reader")?;
    m.add_class::<FchReader>()?;
    m.add_class::<GjfReader>()?;
    m.add_class::<LogReader>()?;
    m.add_class::<MolReader>()?;
    m.add_class::<NonReader>()?;
    m.add_class::<SdfReader>()?;
    m.add_class::<XyzReader>()?;
    m.add_class::<MdeReader>()?;
    m.add_class::<CubReader>()?;
    m.add_class::<DmoReader>()?;
    m.add_class::<InpReader>()?;
    m.add_class::<OutReader>()?;
    m.add_class::<Nics>()?;
    m.add_class::<MfnReader>()?;
    m.add_class::<Marks>()?;
    m.add_class::<Vibe>()?;
    m.add_class::<GauJob>()?;
    m.add_class::<BasData>()?;
    m.add_class::<FileType>()?;
    m.add_function(wrap_pyfunction!(read_lines, &m)?)?;
    m.add_function(wrap_pyfunction!(charge_number, &m)?)?;
    m.add_function(wrap_pyfunction!(get_fileType, &m)?)?;
    m.add_function(wrap_pyfunction!(parse_job, &m)?)?;
    m.add_function(wrap_pyfunction!(shl_basis, &m)?)?;
    m.add_function(wrap_pyfunction!(shl_forms, &m)?)?;
    parent_module.add_submodule(&m)
}
