use crate::base::Mole;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct CubWriter {
    core: rswfn::writer::CubWriter,
}

#[pymethods]
impl CubWriter {
    /// 创建一个新的 CubWriter 实例
    ///
    /// ## 返回
    /// 返回一个 CubWriter 对象
    #[new]
    pub fn new() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::CubWriter::new(),
        };
        Ok(writer)
    }

    /// 设置立方体文件的各个字段
    ///
    /// 所有参数均为可选，仅设置传入的字段。
    ///
    /// ## 参数
    /// - `title`: 标题字符串
    /// - `syms`: 原子符号列表
    /// - `xyzs`: 原子坐标列表
    /// - `obts`: 原子核电荷数列表
    /// - `pos0`: 立方体原点坐标 [x, y, z]
    /// - `size`: 三个方向上的格点数 [nx, ny, nz]
    /// - `step`: 三个方向上的步长 [dx, dy, dz]
    /// - `vals`: 三维格点数据（嵌套 Vec）
    #[pyo3(signature = (title=None, syms=None, xyzs=None, obts=None, pos0=None, size=None, step=None, vals=None))]
    pub fn set(
        &mut self,
        title: Option<String>,
        syms: Option<Vec<String>>,
        xyzs: Option<Vec<[f64; 3]>>,
        obts: Option<Vec<isize>>,
        pos0: Option<[f64; 3]>,
        size: Option<[usize; 3]>,
        step: Option<[f64; 3]>,
        vals: Option<Vec<Vec<f64>>>,
    ) {
        if let Some(title) = title {
            self.core.title = title
        }
        if let Some(syms) = syms {
            self.core.syms = syms
        }
        if let Some(xyzs) = xyzs {
            self.core.xyzs = xyzs;
        }
        if let Some(obts) = obts {
            self.core.obts = obts;
        }
        if let Some(pos0) = pos0 {
            self.core.pos0 = pos0;
        }
        if let Some(size) = size {
            self.core.size = size;
        }
        if let Some(step) = step {
            self.core.step = step;
        }
        if let Some(vals) = vals {
            self.core.vals = vals;
        }
    }

    /// 标题字符串
    #[getter(title)]
    pub fn get_title(&self) -> String {
        self.core.title.clone()
    }

    /// 设置标题字符串
    #[setter(title)]
    pub fn set_title(&mut self, title: String) {
        self.core.title = title;
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

    /// 原子核电荷数列表
    #[getter(obts)]
    pub fn get_obts(&self) -> Vec<isize> {
        self.core.obts.clone()
    }

    /// 设置原子核电荷数列表
    #[setter(obts)]
    pub fn set_obts(&mut self, obts: Vec<isize>) {
        self.core.obts = obts;
    }

    /// 立方体原点坐标
    #[getter(pos0)]
    pub fn get_pos0(&self) -> [f64; 3] {
        self.core.pos0
    }

    /// 设置立方体原点坐标
    #[setter(pos0)]
    pub fn set_pos0(&mut self, pos0: [f64; 3]) {
        self.core.pos0 = pos0;
    }

    /// 三个方向上的格点数
    #[getter(size)]
    pub fn get_size(&self) -> [usize; 3] {
        self.core.size
    }

    /// 设置三个方向上的格点数
    #[setter(size)]
    pub fn set_size(&mut self, size: [usize; 3]) {
        self.core.size = size;
    }

    /// 三个方向上的步长
    #[getter(step)]
    pub fn get_step(&self) -> [f64; 3] {
        self.core.step
    }

    /// 设置三个方向上的步长
    #[setter(step)]
    pub fn set_step(&mut self, step: [f64; 3]) {
        self.core.step = step;
    }

    /// 三维格点数据
    #[getter(vals)]
    pub fn get_vals(&self) -> Vec<Vec<f64>> {
        self.core.vals.clone()
    }

    /// 设置三维格点数据
    #[setter(vals)]
    pub fn set_vals(&mut self, vals: Vec<Vec<f64>>) {
        self.core.vals = vals;
    }

    /// 从 Mole 对象中读取分子信息并填充到 CubWriter
    ///
    /// ## 参数
    /// - `mole`: 包含分子信息的 Mole 对象
    pub fn read_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }

    /// 将立方体数据保存到文件
    ///
    /// ## 参数
    /// - `path`: 保存路径
    pub fn save(&mut self, path: String) {
        self.core.save(&path);
    }
}
