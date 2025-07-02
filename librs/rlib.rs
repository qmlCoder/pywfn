pub mod datas;
pub mod integ;
pub mod march;
pub mod maths;
pub mod matrix;
pub mod props;
pub mod pyfns;
pub mod reads;
pub mod space;
pub mod utils;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;

#[pymodule]
fn rlib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(space::mol_rhos_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::ato_rhos_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::ato_wfns_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::obt_wfns_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::a2m_weits_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::nuc_potential_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::ele_potential_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::get_grids_rs, m)?)?;
    m.add_function(wrap_pyfunction!(space::calc_EIEBA_rs, m)?)?;
    m.add_function(wrap_pyfunction!(march::march_cube_rs, m)?)?;
    m.add_function(wrap_pyfunction!(utils::lag_intpol_rs, m)?)?;
    m.add_function(wrap_pyfunction!(integ::mat_integ_rs, m)?)?;
    m.add_function(wrap_pyfunction!(matrix::ele_mat_rs, m)?)?;

    m.add_function(wrap_pyfunction!(pyfns::search_title, m)?)?;
    m.add_function(wrap_pyfunction!(pyfns::read_geome, m)?)?;
    m.add_function(wrap_pyfunction!(pyfns::read_nbas, m)?)?;
    m.add_function(wrap_pyfunction!(pyfns::read_cmat, m)?)?;
    m.add_function(wrap_pyfunction!(pyfns::read_basis, m)?)?;
    Ok(())
}
