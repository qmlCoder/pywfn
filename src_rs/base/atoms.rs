use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

use crate::base::{Mole, Vec3};

#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct Stm {
    pub core: rswfn::base::Stm,
}

#[pymethods]
impl Stm {
    /// 创建一个新的变换矩阵对象
    ///
    /// **参数**
    /// - `ex`: 矩阵第一行（X轴方向）
    /// - `ey`: 矩阵第二行（Y轴方向）
    /// - `ez`: 矩阵第三行（Z轴方向）
    #[new]
    pub fn new(ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Stm::new(ex, ey, ez),
        })
    }

    /// 设置变换矩阵的三个行向量
    ///
    /// **参数**
    /// - `ex`: 矩阵第一行（X轴方向）
    /// - `ey`: 矩阵第二行（Y轴方向）
    /// - `ez`: 矩阵第三行（Z轴方向）
    pub fn set(&mut self, ex: [f64; 3], ey: [f64; 3], ez: [f64; 3]) {
        // 必须整体重建，否则 mat/non 不会随 ex/ey/ez 一起更新
        self.core = rswfn::base::Stm::new(ex, ey, ez);
    }

    /// 获取矩阵第一行向量
    #[getter]
    pub fn ex(&self) -> [f64; 3] {
        self.core.ex
    }

    /// 获取矩阵第二行向量
    #[getter]
    pub fn ey(&self) -> [f64; 3] {
        self.core.ey
    }

    /// 获取矩阵第三行向量
    #[getter]
    pub fn ez(&self) -> [f64; 3] {
        self.core.ez
    }

    /// 获取完整的 3x3 变换矩阵
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (3, 3)
    #[getter]
    pub fn mat(&self, py: Python) -> Py<PyArray2<f64>> {
        self.core.mat.clone().into_pyarray(py).unbind()
    }

    /// 返回变换矩阵的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}

#[pyclass]
pub struct Atom {
    core: rswfn::base::Atom,
}

#[pymethods]
impl Atom {
    /// 获取原子半径
    #[getter]
    pub fn rad(&self) -> f64 {
        self.core.rad
    }

    /// 获取元素符号
    #[getter]
    pub fn sym(&self) -> String {
        self.core.sym.clone()
    }
}

#[derive(Clone)]
#[pyclass]
pub struct Atoms {
    pub core: rswfn::base::Atoms,
}

#[pymethods]
impl Atoms {
    /// 创建一个新的原子集合
    ///
    /// **参数**
    /// - `atms`: 原子序数列表（例如 6 表示碳）
    /// - `xyzs`: 原子坐标列表，每个坐标为 [x, y, z]
    /// - `crgs`: 原子电荷列表
    #[new]
    pub fn new(atms: Vec<usize>, xyzs: Vec<[f64; 3]>, crgs: Vec<f64>) -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Atoms::new(&atms, &xyzs, &crgs),
        })
    }

    pub fn add(&mut self, atm: usize, xyz: [f64; 3], crg: f64) {
        self.core.add(atm, xyz, crg);
    }

    pub fn get(&self, idx: usize) -> PyResult<Atom> {
        let atom = self.core.get(idx).clone();
        Ok(Atom { core: atom })
    }

    pub fn vec(&self, ai: usize, aj: usize) -> Vec3 {
        Vec3 {
            core: self.core.vec(ai, aj),
        }
    }

    /// 返回原子总数
    pub fn len(&self) -> usize {
        self.core.len()
    }

    pub fn neles(&self) -> [usize; 2] {
        self.core.neles
    }

    pub fn dist(&self, i: usize, j: usize) -> f64 {
        self.core.dist(i, j)
    }

    /// 获取所有原子的元素符号列表
    pub fn syms(&self) -> Vec<String> {
        self.core.syms()
    }

    /// 获取所有原子的坐标列表
    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs().to_vec()
    }

    /// 获取所有原子的原子序数列表
    pub fn nucs(&self) -> Vec<usize> {
        self.core.nucs().clone()
    }

    /// 获取所有原子的电荷列表
    pub fn crgs(&self) -> Vec<f64> {
        self.core.crgs().clone()
    }

    pub fn nebs(&self, atm: usize) -> Vec<usize> {
        self.core.nebs(atm).clone()
    }

    pub fn is_line(&self, atm: usize) -> bool {
        self.core.is_line(atm)
    }

    /// 判断两个原子之间是否存在化学键
    ///
    /// **参数**
    /// - `i`: 第一个原子索引
    /// - `j`: 第二个原子索引
    pub fn is_bond(&self, i: usize, j: usize) -> bool {
        self.core.is_bond(i, j)
    }

    /// 返回原子集合的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }

    pub fn bond_len(&self, ai: usize, aj: usize) -> f64 {
        self.core.bond_len(ai, aj)
    }

    pub fn bond_ang(&self, ai: usize, aj: usize, ak: usize) -> f64 {
        self.core.bond_ang(ai, aj, ak)
    }

    pub fn face_ang(&self, ai: usize, aj: usize, ak: usize, al: usize) -> f64 {
        self.core.face_ang(ai, aj, ak, al)
    }
}
