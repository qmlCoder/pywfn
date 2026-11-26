use crate::base::{Atoms, Basis, Coefs};
use numpy::PyArrayMethods;
use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;
use rswfn;

#[derive(Clone)]
#[pyclass(from_py_object)]
pub struct Mole {
    pub core: rswfn::base::Mole,
}

#[pymethods]
impl Mole {
    /// 从文件加载分子数据
    ///
    /// 支持读取量子化学输出文件（如 .fchk 等格式）并构建分子对象
    ///
    /// **参数**
    /// - `path`: 文件路径字符串
    #[new]
    pub fn new(path: String) -> PyResult<Mole> {
        let mole = Mole {
            core: rswfn::base::Mole::from_file(&path),
        };
        Ok(mole)
    }

    /// 获取分子中的原子集合
    pub fn atoms(&self) -> Atoms {
        let atoms = self.core.atoms().clone();
        Atoms { core: atoms }
    }

    /// 获取分子的基组信息
    pub fn basis(&self) -> Basis {
        let basis = self.core.basis().clone();
        Basis { core: basis }
    }

    /// 获取分子的轨道系数信息
    pub fn coefs(&self) -> Coefs {
        let coefs = self.core.coefs().clone();
        Coefs { core: coefs }
    }

    /// 获取分子的总电荷
    pub fn charge(&self) -> isize {
        self.core.charge()
    }

    /// 获取分子的自旋多重度
    pub fn multip(&self) -> usize {
        self.core.multip()
    }

    /// 获取所有原子的元素符号列表
    pub fn syms(&self) -> Vec<String> {
        self.core.atoms().syms()
    }

    /// 获取所有原子的笛卡尔坐标列表
    pub fn xyzs(&self) -> Vec<[f64; 3]> {
        self.core.atoms().xyzs().clone()
    }

    /// 设置指定基函数形式的轨道系数矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    /// - `cmat`: NumPy 二维数组，轨道系数矩阵
    pub fn set_cmat(&mut self, cmat: Bound<'_, PyArray2<f64>>, form: &str) {
        let cmat = cmat.readonly().as_array().to_owned();
        let mut coefs = self.core.coefs.take().expect("该分子没有系数");
        coefs.set_cmat(form, &cmat, &self.core);
        self.core.coefs = Some(coefs);
    }

    /// 获取指定基函数形式的轨道系数矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_cmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let cmat = self.core.coefs().get_cmat(form, &self.core).clone();
        cmat.into_pyarray(py).unbind()
    }

    /// 获取指定基函数形式的重叠矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_smat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let smat = self.core.basis().get_smat(form, &self.core).clone();
        smat.into_pyarray(py).unbind()
    }

    /// 获取指定基函数形式的密度矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_dmat(&self, py: Python, form: &str) -> Py<PyArray2<f64>> {
        let dmat = self.core.coefs().get_dmat(form, &self.core).clone();
        dmat.into_pyarray(py).unbind()
    }

    /// 判断分子是否为开壳层体系
    pub fn is_open(&self) -> bool {
        self.core.is_open()
    }

    /// 获取分子的边界框
    ///
    /// **返回**
    /// - 包含两个三维坐标的数组：`[min_point, max_point]`
    pub fn border(&self, bord: f64) -> [[f64; 3]; 2] {
        self.core.border(bord)
    }

    /// 返回分子的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }

    pub fn clone(&self) -> Self {
        Mole {
            core: self.core.clone(),
        }
    }
}
