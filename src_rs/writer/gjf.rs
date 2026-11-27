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

    /// 一次性设置 GJF 文件的各个字段
    ///
    /// ## 参数
    /// - `syms`: 原子符号列表
    /// - `xyzs`: 原子坐标列表
    /// - `chk`: checkpoint 文件名
    /// - `job`: 任务类型，如 "opt"、"freq"
    /// - `charge`: 分子总电荷
    /// - `multip`: 自旋多重度
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

    /// checkpoint 文件名
    #[getter(chk)]
    pub fn get_chk(&self) -> String {
        self.core.chk.clone()
    }

    /// 设置 checkpoint 文件名
    #[setter(chk)]
    pub fn set_chk(&mut self, chk: String) {
        self.core.chk = chk;
    }

    /// 任务类型
    #[getter(job)]
    pub fn get_job(&self) -> String {
        self.core.job.clone()
    }

    /// 设置任务类型
    #[setter(job)]
    pub fn set_job(&mut self, job: String) {
        self.core.job = job;
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

    /// 原子符号列表
    #[getter(syms)]
    pub fn get_syms(&self) -> Vec<String> {
        self.core.syms.clone()
    }

    /// 设置原子符号列表
    #[setter(syms)]
    pub fn set_syms(&mut self, syms: Vec<String>) {
        self.core.syms = syms;
    }

    /// 原子坐标列表
    #[getter(xyzs)]
    pub fn get_xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs.clone()
    }

    /// 设置原子坐标列表
    #[setter(xyzs)]
    pub fn set_xyzs(&mut self, xyzs: Vec<[f64; 3]>) {
        self.core.xyzs = xyzs;
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
