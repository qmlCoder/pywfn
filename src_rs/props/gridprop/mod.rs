pub mod density;
pub mod dftgrid;
pub mod potential;
pub mod topolog;
pub mod wfnfunc;

use ndarray::{Array2, Array3};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::prelude::*;
use rswfn;
use rswfn::base::Vec3;

#[pyclass]
pub struct LineGrid {
    pub core: rswfn::props::gridprop::LineGrid,
}

#[pymethods]
impl LineGrid {
    #[new]
    pub fn new() -> Self {
        Self {
            core: rswfn::props::gridprop::LineGrid::new(),
        }
    }

    pub fn set(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64) {
        let p0 = Vec3::from_arr(p0);
        let p1 = Vec3::from_arr(p1);
        self.core.set(&p0, &p1, step);
    }

    pub fn get(&self) -> ([usize; 1], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct RectGrid {
    pub core: rswfn::props::gridprop::RectGrid,
}

#[pymethods]
impl RectGrid {
    #[new]
    pub fn new() -> Self {
        RectGrid {
            core: rswfn::props::gridprop::RectGrid::new(),
        }
    }

    pub fn set_v1(
        &mut self,
        center: [f64; 3],
        normal: [f64; 3],
        vector: [f64; 3],
        size: [f64; 2],
        step: [f64; 2],
    ) {
        let center = Vec3::from_arr(center);
        let normal = Vec3::from_arr(normal);
        let vector = Vec3::from_arr(vector);
        self.core.set_v1(&center, &normal, &vector, size, step);
    }

    pub fn set_v2(&mut self, p0: [f64; 3], p1: [f64; 3], p2: [f64; 3], step: [f64; 2]) {
        let p0 = Vec3::from_arr(p0);
        let p1 = Vec3::from_arr(p1);
        let p2 = Vec3::from_arr(p2);
        self.core.set_v2(&p0, &p1, &p2, step);
    }

    /// 由三个顶点和格点数量确定矩形网格
    ///
    /// ## 参数
    /// - `p00`: 第一个顶点
    /// - `p01`: 第二个顶点
    /// - `p11`: 第三个顶点
    /// - `size`: 两个方向的格点数量
    pub fn set_v3(&mut self, p00: [f64; 3], p01: [f64; 3], p11: [f64; 3], size: [usize; 2]) {
        let p00 = Vec3::from_arr(p00);
        let p01 = Vec3::from_arr(p01);
        let p11 = Vec3::from_arr(p11);
        self.core.set_v3(&p00, &p01, &p11, size);
    }

    /// 计算矩形网格的平面参数
    ///
    /// ## 参数
    /// - `p00`: 第一个顶点
    /// - `p01`: 第二个顶点
    /// - `p11`: 第三个顶点
    ///
    /// ## 返回
    /// 返回 [a, b, c, d] 四个参数
    pub fn paras(&self, p00: [f64; 3], p01: [f64; 3], p11: [f64; 3]) -> [f64; 4] {
        let p00 = Vec3::from_arr(p00);
        let p01 = Vec3::from_arr(p01);
        let p11 = Vec3::from_arr(p11);
        self.core.paras(&p00, &p01, &p11)
    }

    pub fn get(&self) -> ([usize; 2], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct CubeGrid {
    core: rswfn::props::gridprop::CubeGrid,
}

impl CubeGrid {
    /// 由底层对象包装成 Python 对象
    pub fn from_core(core: rswfn::props::gridprop::CubeGrid) -> Self {
        Self { core }
    }
}

#[pymethods]
impl CubeGrid {
    #[new]
    pub fn new() -> Self {
        CubeGrid {
            core: rswfn::props::gridprop::CubeGrid::new(),
        }
    }

    /// 三个方向的格点数
    #[getter]
    pub fn shape(&self) -> [usize; 3] {
        self.core.shape
    }

    /// 格点起点坐标
    #[getter]
    pub fn start(&self) -> [f64; 3] {
        self.core.start
    }

    /// 三个方向的格点步长
    #[getter]
    pub fn steps(&self) -> [f64; 3] {
        self.core.steps
    }

    pub fn set_v1(&mut self, p0: [f64; 3], p1: [f64; 3], step: f64, bord: f64) {
        self.core.set_v1(&p0, &p1, step, bord);
    }

    pub fn get(&self, py: Python) -> ([usize; 3], Py<PyArray2<f64>>) {
        let (shape, grids) = self.core.get();
        let arr = Array2::from_shape_vec((grids.len(), 3), grids.concat()).unwrap();
        (shape, arr.into_pyarray(py).unbind())
    }
}

/// 格点上函数值及其各级导数
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct GridVals {
    pub core: rswfn::props::gridprop::GridVals,
}

#[pymethods]
impl GridVals {
    /// 创建一个格点函数值对象
    ///
    /// ## 参数
    /// - `val0`: 函数值
    /// - `val1`: 梯度 [dx, dy, dz]
    /// - `val2`: Hessian 矩阵 3x3
    #[new]
    #[pyo3(signature = (val0=0.0, val1=None, val2=None))]
    pub fn new(val0: f64, val1: Option<[f64; 3]>, val2: Option<[[f64; 3]; 3]>) -> Self {
        let val1 = val1.unwrap_or([0.0, 0.0, 0.0]);
        let val2 = val2.unwrap_or([[0.0; 3]; 3]);
        Self {
            core: rswfn::props::gridprop::GridVals::new(val0, val1, val2),
        }
    }

    /// 设置函数值及其各级导数
    ///
    /// ## 参数
    /// - `val0`: 函数值
    /// - `val1`: 梯度 [dx, dy, dz]
    /// - `val2`: Hessian 矩阵 3x3
    pub fn set(&mut self, val0: f64, val1: [f64; 3], val2: [[f64; 3]; 3]) {
        self.core.set(val0, val1, val2);
    }

    /// 返回一个全零的格点函数值对象
    #[staticmethod]
    pub fn non() -> Self {
        Self {
            core: rswfn::props::gridprop::GridVals::non(),
        }
    }

    /// 函数值
    #[getter]
    pub fn val0(&self) -> f64 {
        self.core.val0
    }

    /// 梯度
    #[getter]
    pub fn val1(&self) -> [f64; 3] {
        self.core.val1
    }

    /// Hessian 矩阵
    #[getter]
    pub fn val2(&self) -> [[f64; 3]; 3] {
        self.core.val2
    }

    /// 梯度的模
    pub fn grad_norm(&self) -> f64 {
        self.core.grad_norm()
    }

    /// 梯度模的平方
    pub fn grad_norm2(&self) -> f64 {
        self.core.grad_norm2()
    }

    /// Hessian 矩阵展平为长度 9 的列表
    pub fn val2s(&self) -> [f64; 9] {
        self.core.val2s()
    }
}

/// 把 rswfn 的格点函数值包装成 Python 对象
pub fn wrap_gridvals(core: rswfn::props::gridprop::GridVals) -> GridVals {
    GridVals { core }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "gridprop")?;
    density::register_module(&m)?;
    dftgrid::register_module(&m)?;
    wfnfunc::register_module(&m)?;
    potential::register_module(&m)?;
    topolog::register_module(&m)?;
    m.add_class::<LineGrid>()?;
    m.add_class::<RectGrid>()?;
    m.add_class::<CubeGrid>()?;
    m.add_class::<GridVals>()?;
    parent_module.add_submodule(&m)
}

/// 将 `GridVals` 列表按值、梯度、Hessian 拆分为三个独立的向量。
///
/// * `vals` - 网格值列表
/// * `level` - 计算级别（0=仅值，1=值+梯度，2=值+梯度+Hessian）
///
/// 返回 `(values, gradients, hessians)` 三元组。
fn vals_split(
    vals: &Vec<rswfn::props::gridprop::GridVals>,
    level: usize,
) -> (Vec<f64>, Vec<Grad>, Vec<Hess>) {
    let mut vals0 = vec![];
    let mut vals1 = vec![];
    let mut vals2 = vec![];
    for val in vals {
        vals0.push(val.val0);
        if level > 0 {
            vals1.push(val.val1);
        }
        if level > 1 {
            vals2.push(val.val2);
        }
    }
    (vals0, vals1, vals2)
}

type Grad = [f64; 3];
type Hess = [[f64; 3]; 3];
/// 将梯度向量列表转换为形状为 `(n, 3)` 的二维数组。
pub fn grads2array(grads: Vec<Grad>) -> Array2<f64> {
    let ngrad = grads.len();
    let datas = grads.into_flattened();
    Array2::from_shape_vec((ngrad, 3), datas).unwrap()
}
/// 将 Hessian 矩阵列表转换为形状为 `(n, 3, 3)` 的三维数组。
pub fn hesss2array(hesss: Vec<Hess>) -> Array3<f64> {
    let nhess = hesss.len();
    let datas = hesss.into_flattened().into_flattened();
    Array3::from_shape_vec((nhess, 3, 3), datas).unwrap()
}

pub fn vals2numpy(
    py: Python,
    vals: &Vec<rswfn::props::gridprop::GridVals>,
    level: usize,
) -> (Py<PyArray1<f64>>, Py<PyArray2<f64>>, Py<PyArray3<f64>>) {
    let (vals, grad, hess) = vals_split(vals, level);
    let vals = vals.into_pyarray(py).unbind();
    let grad = grads2array(grad).into_pyarray(py).unbind();
    let hess = hesss2array(hess).into_pyarray(py).unbind();
    (vals, grad, hess)
}
