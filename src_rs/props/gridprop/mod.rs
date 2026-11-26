pub mod density;
pub mod dftgrid;
pub mod potential;
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

    pub fn get(&self) -> ([usize; 2], Vec<[f64; 3]>) {
        let (shape, grids) = self.core.get();
        (shape, grids.clone())
    }
}

#[pyclass]
pub struct CubeGrid {
    core: rswfn::props::gridprop::CubeGrid,
}

#[pymethods]
impl CubeGrid {
    #[new]
    pub fn new() -> Self {
        CubeGrid {
            core: rswfn::props::gridprop::CubeGrid::new(),
        }
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

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "gridprop")?;
    density::register_module(&m)?;
    dftgrid::register_module(&m)?;
    wfnfunc::register_module(&m)?;
    potential::register_module(&m)?;
    m.add_class::<LineGrid>()?;
    m.add_class::<RectGrid>()?;
    m.add_class::<CubeGrid>()?;
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
