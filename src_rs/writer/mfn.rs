use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct MfnWriter {
    core: rswfn::writer::MfnWriter,
}

#[pymethods]
impl MfnWriter {
    /// 创建一个新的 MfnWriter 实例
    ///
    /// ## 返回
    /// 返回一个 MfnWriter 对象
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::MfnWriter::new(),
        };
        Ok(writer)
    }

    /// 从 Mole 对象构建 MfnWriter
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    ///
    /// ## 返回
    /// 返回填好分子信息的 MfnWriter 对象
    #[staticmethod]
    pub fn from_mole(mole: Mole) -> Self {
        Self {
            core: rswfn::writer::MfnWriter::from_mole(&mole.core),
        }
    }

    /// 分子对象
    #[getter(mole)]
    pub fn get_mole(&self) -> Mole {
        Mole {
            core: self.core.mole.clone(),
        }
    }

    /// 设置分子对象
    #[setter(mole)]
    pub fn set_mole(&mut self, mole: PyRef<'_, Mole>) {
        self.core.mole = mole.core.clone();
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

    /// 生成 mfn 文件的文本内容
    ///
    /// ## 返回
    /// 返回格式化的 mfn 文件字符串
    pub fn build(&self) -> String {
        self.core.build()
    }

    /// 将 mfn 文件内容保存到磁盘
    ///
    /// ## 参数
    /// - `path`: 保存路径
    pub fn save(&self, path: String) {
        self.core.save(&path);
    }

    /// 从分子对象中读取数据
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn read_mole(&mut self, mole: Mole) {
        self.core.mole = mole.core;
    }
}
