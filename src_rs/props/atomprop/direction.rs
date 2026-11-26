use pyo3::prelude::*;

use rswfn;

use crate::base::{Mole, Stm};

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::props::atomprop::direction::Calculator<'_> {
        rswfn::props::atomprop::direction::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的原子方向计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算原子的法向量
    ///
    /// 返回指定原子的法向量方向，用于描述原子在分子中的空间取向。
    ///
    /// ## 参数
    /// - `atm`: 原子索引
    ///
    /// ## 返回
    /// 法向量 [x, y, z]，如果无法计算则返回 None
    pub fn normal_vector(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().normal_vector(atm);
        match dir {
            Some(dir) => Some(dir.into_arr()),
            None => None,
        }
    }

    /// 获取原子的归一化法向量
    ///
    /// 与 `normal_vector` 功能类似，计算原子的归一化法向量。
    ///
    /// ## 参数
    /// - `atm`: 原子索引
    ///
    /// ## 返回
    /// 归一化法向量 [x, y, z]，如果无法计算则返回 None
    pub fn get_normal(&self, atm: usize) -> Option<[f64; 3]> {
        let dir = self.caler().get_normal(atm);
        match dir {
            Some(dir) => Some(dir.into_arr()),
            None => None,
        }
    }

    /// 计算局域坐标系（Local Coordinate System）
    ///
    /// 基于指定原子及其邻域原子构建局域坐标系旋转矩阵。
    ///
    /// ## 参数
    /// - `atm`: 中心原子索引
    /// - `neb`: 可选的邻域原子索引，用于确定坐标系方向
    ///
    /// ## 返回
    /// 局域坐标系旋转矩阵（Stm），如果无法构建则返回 None
    pub fn LCS(&mut self, atm: usize, neb: Option<usize>) -> Option<Stm> {
        let stm = self.caler().LCS(atm, neb);
        match stm {
            Some(stm) => Some(Stm { core: stm }),
            None => None,
        }
    }

    /// 计算原子的反应方向
    ///
    /// 返回指定原子可能发生化学反应的方向向量，用于预测反应位点。
    ///
    /// ## 参数
    /// - `atm`: 原子索引
    ///
    /// ## 返回
    /// 多个反应方向向量 [x, y, z] 的列表
    pub fn reactions(&self, atm: usize) -> Vec<[f64; 3]> {
        let dirs = self.caler().reactions(atm);
        let dirs: Vec<[f64; 3]> = dirs.iter().map(|dir| dir.into_arr()).collect();
        dirs
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "direction")?;
    m.add_class::<Calculator>()?;
    parent_module.add_submodule(&m)
}

#[pymodule]
pub mod direction {
    #[pymodule_export]
    use super::Calculator;
}
