use crate::base::{Atoms, Basis, Coefs, Mole};
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct FchWriter {
    core: rswfn::writer::FchWriter,
}

#[pymethods]
impl FchWriter {
    /// 创建一个新的 FchWriter 实例
    ///
    /// ## 返回
    /// 返回一个 FchWriter 对象
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::FchWriter::new(),
        };
        Ok(writer)
    }

    /// 从 Mole 对象构建 FchWriter
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 返回填好分子信息的 FchWriter 对象
    #[staticmethod]
    pub fn from_mole(mole: Mole) -> Self {
        Self {
            core: rswfn::writer::FchWriter::from_mole(&mole.core),
        }
    }

    /// 分子总电荷
    #[getter(charge)]
    pub fn get_charge(&self) -> isize {
        self.core.charge
    }

    /// 设置分子总电荷
    #[setter(charge)]
    pub fn set_charge(&mut self, charge: isize) {
        self.core.charge = charge;
    }

    /// 自旋多重度
    #[getter(multip)]
    pub fn get_multip(&self) -> usize {
        self.core.multip
    }

    /// 设置自旋多重度
    #[setter(multip)]
    pub fn set_multip(&mut self, multip: usize) {
        self.core.multip = multip;
    }

    /// 分子能量
    #[getter(energy)]
    pub fn get_energy(&self) -> f64 {
        self.core.energy
    }

    /// 设置分子能量
    #[setter(energy)]
    pub fn set_energy(&mut self, energy: f64) {
        self.core.energy = energy;
    }

    /// 原子信息
    #[getter(atoms)]
    pub fn get_atoms(&self) -> Atoms {
        Atoms {
            core: self.core.atoms.clone(),
        }
    }

    /// 设置原子信息
    #[setter(atoms)]
    pub fn set_atoms(&mut self, atoms: PyRef<'_, Atoms>) {
        self.core.atoms = atoms.core.clone();
    }

    /// 基组信息
    #[getter(basis)]
    pub fn get_basis(&self) -> Basis {
        Basis {
            core: self.core.basis.clone(),
        }
    }

    /// 设置基组信息
    #[setter(basis)]
    pub fn set_basis(&mut self, basis: PyRef<'_, Basis>) {
        self.core.basis = basis.core.clone();
    }

    /// 轨道系数信息
    #[getter(coefs)]
    pub fn get_coefs(&self) -> Coefs {
        Coefs {
            core: self.core.coefs.clone(),
        }
    }

    /// 设置轨道系数信息
    #[setter(coefs)]
    pub fn set_coefs(&mut self, coefs: PyRef<'_, Coefs>) {
        self.core.coefs = coefs.core.clone();
    }

    /// 生成 fchk 文件的文本内容
    ///
    /// ## 返回
    /// 返回格式化的 fchk 文件字符串
    pub fn build(&self) -> String {
        self.core.build()
    }

    /// 将 fchk 文件内容保存到磁盘
    ///
    /// ## 参数
    /// - `path`: 保存路径
    pub fn save(&self, path: String) {
        self.core.save(&path);
    }

    /// 从 Mole 对象中读取分子信息
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn read_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }
}
