pub mod march;
pub mod datas;
pub mod space;
pub mod utils;
pub mod integ;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pymodule]
fn rlib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(space::mol_rhos, m)?)?;
    m.add_function(wrap_pyfunction!(space::obt_wfns, m)?)?;
    m.add_function(wrap_pyfunction!(space::a2m_weits, m)?)?;
    m.add_function(wrap_pyfunction!(march::march_cube, m)?)?;
    m.add_function(wrap_pyfunction!(utils::lag_intpol_rs, m)?)?;
    m.add_function(wrap_pyfunction!(integ::mat_integ_rs, m)?)?;
    Ok(())
}
