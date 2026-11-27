use pyo3::prelude::*;
use rswfn;

/// 玻尔半径（单位：埃）
#[pyfunction]
fn BOHR() -> f64 {
    rswfn::datas::consts::BOHR
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "consts")?;
    m.add_function(wrap_pyfunction!(BOHR, &m)?)?;
    parent_module.add_submodule(&m)
}
