#![allow(non_snake_case)]

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;
use rswfn;

/// 一次函数求值的结果
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct FuncEval {
    pub core: rswfn::maths::optim::FuncEval,
}

#[pymethods]
impl FuncEval {
    /// 创建一次函数求值结果
    ///
    /// ## 参数
    /// - `value`: 函数值
    /// - `grad`: 梯度列表，可选
    /// - `hessian`: Hessian 矩阵，可选
    #[new]
    #[pyo3(signature = (value, grad=None, hessian=None))]
    pub fn new(
        value: f64,
        grad: Option<Vec<f64>>,
        hessian: Option<Bound<'_, PyArray2<f64>>>,
    ) -> Self {
        let hessian = hessian.map(|h| h.readonly().as_array().to_owned());
        Self {
            core: rswfn::maths::optim::FuncEval::new(value, grad, hessian),
        }
    }

    /// 仅给定函数值
    #[staticmethod]
    pub fn from_value(value: f64) -> Self {
        Self {
            core: rswfn::maths::optim::FuncEval::from_value(value),
        }
    }

    /// 给定函数值和梯度
    #[staticmethod]
    pub fn from_grad(value: f64, grad: Vec<f64>) -> Self {
        Self {
            core: rswfn::maths::optim::FuncEval::from_grad(value, grad),
        }
    }

    /// 给定函数值和 Hessian 矩阵
    #[staticmethod]
    pub fn from_hessian(value: f64, hessian: Bound<'_, PyArray2<f64>>) -> Self {
        let hessian = hessian.readonly().as_array().to_owned();
        Self {
            core: rswfn::maths::optim::FuncEval::from_hessian(value, hessian),
        }
    }

    /// 函数值
    #[getter]
    pub fn value(&self) -> f64 {
        self.core.value
    }

    /// 梯度
    #[getter]
    pub fn grad(&self) -> Option<Vec<f64>> {
        self.core.grad.clone()
    }

    /// Hessian 矩阵
    #[getter]
    pub fn hessian(&self, py: Python) -> Option<Py<PyArray2<f64>>> {
        self.core
            .hessian
            .clone()
            .map(|hess| hess.into_pyarray(py).unbind())
    }
}

/// 把 Python 返回值转换为一次函数求值结果
///
/// 支持三种写法：直接返回浮点数、返回 `FuncEval`、返回 `(value, grad, hessian)` 元组。
fn parse_funceval(py: Python, obj: &Bound<'_, PyAny>) -> rswfn::maths::optim::FuncEval {
    if let Ok(value) = obj.extract::<f64>() {
        return rswfn::maths::optim::FuncEval::from_value(value);
    }
    if let Ok(feval) = obj.extract::<FuncEval>() {
        return feval.core;
    }
    let tup = obj
        .cast::<PyTuple>()
        .expect("优化回调应返回浮点数、FuncEval 或 (value, grad, hessian) 元组");
    let value: f64 = tup.get_item(0).unwrap().extract().unwrap();
    let grad: Option<Vec<f64>> = match tup.len() {
        0 | 1 => None,
        _ => tup
            .get_item(1)
            .ok()
            .and_then(|item| item.extract::<Vec<f64>>().ok()),
    };
    let hessian: Option<Array2<f64>> = match tup.len() {
        0..=2 => None,
        _ => tup
            .get_item(2)
            .ok()
            .and_then(|item| item.extract::<Bound<'_, PyArray2<f64>>>().ok())
            .map(|arr| arr.readonly().as_array().to_owned()),
    };
    let _ = py;
    rswfn::maths::optim::FuncEval::new(value, grad, hessian)
}

macro_rules! optim_mod {
    ($name:ident) => {
        pub mod $name {
            use super::*;

            /// 最小化目标函数
            ///
            /// ## 参数
            /// - `func`: Python 回调，接收参数列表，返回函数值（或 `FuncEval`、`(value, grad, hessian)`）
            /// - `paras`: 初始参数列表
            /// - `conv_thresh`: 收敛阈值
            ///
            /// ## 返回
            /// 返回优化后的参数列表
            #[pyfunction]
            pub fn minimize(
                py: Python,
                func: Py<PyAny>,
                paras: Vec<f64>,
                conv_thresh: f64,
            ) -> PyResult<Vec<f64>> {
                let callback = func.clone_ref(py);
                let res = py.detach(move || {
                    rswfn::maths::optim::$name::minimize(
                        |p: &Vec<f64>| {
                            Python::attach(|py| {
                                let args = (p.clone(),);
                                let out = callback
                                    .call1(py, args)
                                    .expect("优化回调执行失败");
                                parse_funceval(py, out.bind(py))
                            })
                        },
                        &paras,
                        conv_thresh,
                    )
                });
                res.map_err(PyValueError::new_err)
            }

            pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
                let m = PyModule::new(parent_module.py(), stringify!($name))?;
                m.add_function(wrap_pyfunction!(minimize, &m)?)?;
                parent_module.add_submodule(&m)
            }
        }
    };
}

optim_mod!(cg);
optim_mod!(dn);
optim_mod!(nt);
optim_mod!(qn);
optim_mod!(sd);

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "optim")?;
    m.add_class::<FuncEval>()?;
    cg::register_module(&m)?;
    dn::register_module(&m)?;
    nt::register_module(&m)?;
    qn::register_module(&m)?;
    sd::register_module(&m)?;
    parent_module.add_submodule(&m)
}
