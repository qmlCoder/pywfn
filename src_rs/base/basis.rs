use numpy::{IntoPyArray, PyArray2};
use pyo3::prelude::*;

use rswfn;

use crate::base::Mole;

#[pyclass]
pub struct Atos {
    pub core: rswfn::base::Atos,
}

#[pymethods]
impl Atos {
    /// 返回指定基函数形式的原子轨道数量
    ///
    /// **参数**
    /// - `form`: 基函数形式，如 "STO-3G"
    pub fn len(&self, form: &str) -> usize {
        self.core.len(form)
    }

    /// 获取指定基函数形式下各壳层对应的原子索引列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn atms(&self, form: &str) -> Vec<usize> {
        self.core.atms(form)
    }

    /// 获取指定基函数形式下各壳层对应的壳层索引列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn shls(&self, form: &str) -> Vec<usize> {
        self.core.shls(form)
    }

    /// 获取指定基函数形式下各壳层的符号表示列表
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn syms(&self, form: &str) -> Vec<String> {
        self.core.syms(form)
    }
}

#[pyclass]
pub struct Shl {
    pub core: rswfn::base::Shl,
}

#[pymethods]
impl Shl {
    /// 获取该壳层所属的原子索引
    #[getter]
    pub fn atm(&self) -> usize {
        self.core.atm
    }
    /// 获取该壳层的索引编号
    #[getter]
    pub fn shl(&self) -> usize {
        self.core.shl
    }
    /// 获取该壳层的角动量量子数
    #[getter]
    pub fn ang(&self) -> usize {
        self.core.ang
    }
    /// 获取该壳层的起始轨道索引位置
    #[getter]
    pub fn loc(&self) -> usize {
        self.core.loc()
    }
    /// 获取该壳层的基函数形式
    #[getter]
    pub fn form(&self) -> String {
        self.core.form.clone()
    }
    /// 获取该壳层的高斯指数（alpha 值）列表
    #[getter]
    pub fn alps(&self) -> Vec<f64> {
        self.core.alps.clone()
    }
    /// 获取该壳层的收缩系数列表
    #[getter]
    pub fn coes(&self) -> Vec<f64> {
        self.core.coes.clone()
    }
}

#[pyclass]
pub struct Shls {
    pub core: rswfn::base::Shls,
}

#[pymethods]
impl Shls {
    /// 获取指定索引的壳层
    ///
    /// **参数**
    /// - `i`: 壳层索引
    ///
    /// **返回**
    /// - 壳层对象
    pub fn get(&self, i: usize) -> Shl {
        Shl {
            core: self.core.get(i).clone(),
        }
    }

    /// 返回壳层总数
    pub fn len(&self) -> usize {
        self.core.len()
    }

    /// 获取指定基函数形式的轨道上三角索引对列表
    ///
    /// 返回每对轨道的 (i, j) 索引，用于构建重叠矩阵等
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn uplo(&self, form: &str) -> Vec<[usize; 2]> {
        self.core.uplo(form).clone()
    }
}

#[derive(Clone)]
#[pyclass]
pub struct Basis {
    pub core: rswfn::base::Basis,
}

#[pymethods]
impl Basis {
    /// 返回基组的字符串表示
    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }

    /// 获取原子轨道壳层索引信息
    pub fn atos(&self) -> Atos {
        Atos {
            core: self.core.atos().clone(),
        }
    }

    /// 获取所有壳层对象
    pub fn shls(&self) -> Shls {
        Shls {
            core: self.core.shls.clone(),
        }
    }

    /// 获取指定基函数形式的重叠矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式
    /// - `mole`: 分子对象
    ///
    /// **返回**
    /// - NumPy 二维数组，形状为 (n_orb, n_orb)
    pub fn get_smat(&self, py: Python, form: &str, mole: Mole) -> Py<PyArray2<f64>> {
        let smat = self.core.get_smat(form, &mole.core).to_owned();
        smat.into_pyarray(py).unbind()
    }

    /// 获取指定基函数形式的原子对索引列表
    ///
    /// 返回每对相互作用的原子索引 (i, j)，用于构建原子基的重叠矩阵
    ///
    /// **参数**
    /// - `form`: 基函数形式
    pub fn atm_uplo(&self, form: &str) -> Vec<[usize; 2]> {
        let uplo = self.core.atm_uplo(form).to_owned();
        uplo
    }
}
