use pyo3::prelude::*;

use rswfn;

#[pyclass]
struct Tool {
    inner: rswfn::tools::log::Tool,
}

#[pymethods]
impl Tool {
    #[new]
    fn new() -> Tool {
        Tool {
            inner: rswfn::tools::log::Tool::new(),
        }
    }

    pub fn get_SI(&self, paths: Vec<String>) -> String {
        self.inner.get_SI(&paths)
    }

    pub fn split_opt(&self, path: String) -> Vec<(Vec<u32>, Vec<[f64; 3]>, f64)> {
        self.inner.split_opt(path)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "log")?; // 为当前的rs文件创建一个子模块
    m.add_class::<Tool>()?; // 添加这个rs文件中的函数到子模块中
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
