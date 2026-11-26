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
