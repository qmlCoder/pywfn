use pyo3::prelude::*;

use numpy::PyReadonlyArray2;
use rswfn;

use crate::utils;

#[pyclass]
pub struct Coefs {
    pub core: rswfn::base::Coefs,
}

#[pymethods]
impl Coefs {
    #[new]
    fn new(
        ato_atms: Vec<usize>,
        ato_shls: Vec<usize>,
        ato_syms: Vec<String>,
        obt_engs: Vec<f64>,
        obt_occs: Vec<usize>,
        cmat: PyReadonlyArray2<f64>, // 传入numpy矩阵
    ) -> Coefs {
        Coefs {
            core: rswfn::base::Coefs::new(
                ato_atms,
                ato_shls,
                ato_syms,
                obt_engs,
                obt_occs,
                cmat.as_array().to_owned(),
            ),
        }
    }

    pub fn get_cmat(&self, form: String) -> Vec<Vec<f64>> {
        let cmat = self.core.get_cmat(&form);
        utils::arr2vec(cmat)
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.core)
    }
}
