use pyo3::prelude::*;
use rswfn;

#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct Vec3 {
    pub core: rswfn::base::Vec3,
}

#[pymethods]
impl Vec3 {
    /// 创建一个新的三维向量
    #[new]
    pub fn new(x: f64, y: f64, z: f64) -> PyResult<Self> {
        Ok(Self {
            core: rswfn::base::Vec3::new(x, y, z),
        })
    }

    /// 从长度为3的数组创建三维向量
    #[staticmethod]
    pub fn from_arr(vec: [f64; 3]) -> Self {
        Self {
            core: rswfn::base::Vec3::from_arr(vec),
        }
    }

    /// 从列表中创建三维向量
    #[staticmethod]
    pub fn from_vec(vec: Vec<f64>) -> Self {
        Self {
            core: rswfn::base::Vec3::from_vec(&vec),
        }
    }

    /// 返回零向量 (0, 0, 0)
    #[staticmethod]
    pub fn zero() -> Self {
        Self {
            core: rswfn::base::Vec3::zero(),
        }
    }

    /// 生成一个随机单位向量
    #[staticmethod]
    pub fn rand() -> Self {
        Self {
            core: rswfn::base::Vec3::rand(),
        }
    }

    /// 将向量转换为长度为3的数组
    pub fn into_arr(&self) -> [f64; 3] {
        self.core.into_arr()
    }

    /// 将向量转换为列表
    pub fn into_vec(&self) -> Vec<f64> {
        self.core.into_vec()
    }

    /// x 分量
    #[getter]
    pub fn x(&self) -> f64 {
        self.core.into_arr()[0]
    }

    /// y 分量
    #[getter]
    pub fn y(&self) -> f64 {
        self.core.into_arr()[1]
    }

    /// z 分量
    #[getter]
    pub fn z(&self) -> f64 {
        self.core.into_arr()[2]
    }

    /// 返回归一化后的单位向量
    pub fn normalize(&self) -> Vec3 {
        Self {
            core: self.core.normalize(),
        }
    }

    /// 返回向量的长度（模）
    pub fn length(&self) -> f64 {
        self.core.length()
    }

    /// 计算与另一个向量的叉积
    pub fn cross(&self, other: Vec3) -> Vec3 {
        Self {
            core: self.core.cross(&other.core),
        }
    }

    /// 计算与另一个向量的点积
    pub fn dot(&self, other: Vec3) -> f64 {
        self.core.dot(&other.core)
    }

    /// 两向量之间的夹角，范围 0-PI
    pub fn angle(&self, other: Vec3) -> f64 {
        self.core.angle(&other.core)
    }

    /// 计算向量在另一个向量方向上的投影
    pub fn project_onto(&self, other: Vec3) -> Vec3 {
        Self {
            core: self.core.project_onto(&other.core),
        }
    }

    /// 绕指定轴旋转一个点，self 代表 point
    pub fn rotate_arrow(&self, axis: Vec3, angle: f64) -> Vec3 {
        Self {
            core: self.core.rotate_arrow(&axis.core, angle),
        }
    }

    /// 绕指定中心点和轴旋转一个点，self 代表 point
    pub fn rotate_point(&self, cent: Vec3, axis: Vec3, angle: f64) -> Vec3 {
        Self {
            core: self.core.rotate_point(&cent.core, &axis.core, angle),
        }
    }

    /// 计算点到平面的距离，平面由 ax + by + cz + d = 0 定义
    pub fn dist_to_plan(&self, plan: [f64; 4]) -> f64 {
        self.core.dist_to_plan(&plan)
    }

    pub fn __add__(&self, other: Vec3) -> Vec3 {
        Self {
            core: &self.core + &other.core,
        }
    }

    pub fn __radd__(&self, other: Vec3) -> Vec3 {
        Self {
            core: &other.core + &self.core,
        }
    }

    pub fn __sub__(&self, other: Vec3) -> Vec3 {
        Self {
            core: &self.core - &other.core,
        }
    }

    pub fn __rsub__(&self, other: Vec3) -> Vec3 {
        Self {
            core: &other.core - &self.core,
        }
    }

    pub fn __mul__(&self, scalar: f64) -> Vec3 {
        Self {
            core: &self.core * scalar,
        }
    }

    pub fn __rmul__(&self, scalar: f64) -> Vec3 {
        Self {
            core: &self.core * scalar,
        }
    }

    pub fn __truediv__(&self, scalar: f64) -> Vec3 {
        Self {
            core: self.core.clone() / scalar,
        }
    }

    pub fn __neg__(&self) -> Vec3 {
        Self {
            core: &self.core * -1.0,
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
