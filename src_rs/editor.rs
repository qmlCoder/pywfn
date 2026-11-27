#![allow(non_snake_case)]

use pyo3::prelude::*;
use rswfn;

use crate::base::Mole;

/// 旋转所有给定点
///
/// ## 参数
/// - `points`: 待旋转的点列表
/// - `cent`: 旋转中心
/// - `axis`: 旋转轴方向
/// - `angle`: 旋转角度（弧度）
///
/// ## 返回
/// - 旋转后的点列表
#[pyfunction]
pub fn rotate_points(
    points: Vec<[f64; 3]>,
    cent: [f64; 3],
    axis: [f64; 3],
    angle: f64,
) -> Vec<[f64; 3]> {
    rswfn::editor::rotate_points(&points, &cent, &axis, angle)
}

#[pyclass]
#[derive(Clone)]
pub struct Editor {
    pub mole: Mole,
}

#[pymethods]
impl Editor {
    /// 创建一个分子编辑器
    ///
    /// ## 参数
    /// - `mole`: 需要编辑的分子对象
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 获取与指定原子直接相连的原子分组
    ///
    /// ## 参数
    /// - `atm`: 目标原子索引
    /// - `excludes`: 需要排除的原子索引
    ///
    /// ## 返回
    /// - 相连原子的索引列表
    pub fn get_group(&self, atm: usize, excludes: Vec<usize>) -> Vec<usize> {
        let editor = rswfn::editor::Editor::new(&self.mole.core);
        editor.get_group(atm, &excludes)
    }

    /// 按键长设置分子结构
    ///
    /// ## 参数
    /// - `a1`: 第一个原子索引
    /// - `a2`: 第二个原子索引
    /// - `value`: 目标键长
    /// - `stype`: 结构类型
    ///
    /// ## 返回
    /// - 修改后的新分子对象
    pub fn set_bl_by(&self, a1: usize, a2: usize, value: f64, stype: usize) -> Mole {
        let editor = rswfn::editor::Editor::new(&self.mole.core);
        Mole {
            core: editor.set_bl_by(a1, a2, value, stype),
        }
    }

    /// 按键角设置分子结构
    ///
    /// ## 参数
    /// - `a1`: 第一个原子索引
    /// - `a2`: 第二个原子索引
    /// - `a3`: 第三个原子索引
    /// - `value`: 目标键角（弧度）
    /// - `btype`: 结构类型
    ///
    /// ## 返回
    /// - 修改后的新分子对象
    pub fn set_ba_by(&self, a1: usize, a2: usize, a3: usize, value: f64, btype: usize) -> Mole {
        let editor = rswfn::editor::Editor::new(&self.mole.core);
        Mole {
            core: editor.set_ba_by(a1, a2, a3, value, btype),
        }
    }

    /// 按二面角设置分子结构
    ///
    /// ## 参数
    /// - `a1`: 第一个原子索引
    /// - `a2`: 第二个原子索引
    /// - `value`: 目标二面角（弧度）
    /// - `rtype`: 结构类型
    ///
    /// ## 返回
    /// - 修改后的新分子对象
    pub fn set_da_by(&self, a1: usize, a2: usize, value: f64, rtype: usize) -> Mole {
        let editor = rswfn::editor::Editor::new(&self.mole.core);
        Mole {
            core: editor.set_da_by(a1, a2, value, rtype),
        }
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "editor")?;
    m.add_class::<Editor>()?;
    m.add_function(wrap_pyfunction!(rotate_points, &m)?)?;
    parent_module.add_submodule(&m)
}
