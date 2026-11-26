use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct GjfWriter {
    core: rswfn::writer::GjfWriter,
}

#[pymethods]
impl GjfWriter {
    /// 创建一个新的 GjfWriter 实例
    ///
    /// ## 返回
    /// 返回一个 GjfWriter 对象
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::GjfWriter::new(),
        };
        Ok(writer)
    }

    pub fn set(
        &mut self,
        syms: Vec<String>,
        xyzs: Vec<[f64; 3]>,
        chk: String,
        job: String,
        charge: isize,
        multip: usize,
    ) {
        self.core.syms = syms;
        self.core.xyzs = xyzs;
        self.core.chk = chk;
        self.core.job = job;
        self.core.charge = charge;
        self.core.multip = multip;
    }

    /// 生成 GJF 文件的文本内容
    ///
    /// ## 返回
    /// 返回格式化的 GJF 文件字符串
    pub fn build(&self) -> String {
        self.core.build()
    }

    /// 将 GJF 文件内容保存到磁盘
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
        self.core.read_mole(&mole.core);
    }
}
