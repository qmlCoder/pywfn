use crate::base::Mole;
use ndarray::Array2;
use pyo3::prelude::*;
use rswfn;
use rswfn::writer::Writer;

#[pyclass]
pub struct DmoWriter {
    core: rswfn::writer::DmoWriter,
}

/// 把二维 `Vec<Vec<f64>>` 转成 `Array2`，行长度不一致时返回 Python 的 ValueError。
fn to_array2(rows: &[Vec<f64>], name: &str) -> PyResult<Array2<f64>> {
    let nrows = rows.len();
    let ncols = rows.first().map_or(0, |row| row.len());
    if rows.iter().any(|row| row.len() != ncols) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "{name} 的每一行长度必须一致"
        )));
    }
    let flat: Vec<f64> = rows.iter().flatten().copied().collect();
    Ok(Array2::from_shape_vec((nrows, ncols), flat)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?)
}

#[pymethods]
impl DmoWriter {
    #[new]
    pub fn blank() -> PyResult<Self> {
        let writer = Self {
            core: rswfn::writer::DmoWriter::blank(),
        };
        Ok(writer)
    }

    #[getter(syms)]
    pub fn get_syms(&self) -> Vec<String> {
        self.core.syms.clone()
    }

    #[setter(syms)]
    pub fn set_syms(&mut self, syms: Vec<String>) {
        self.core.syms = syms;
    }

    #[getter(xyzs)]
    pub fn get_xyzs(&self) -> Vec<[f64; 3]> {
        self.core.xyzs.clone()
    }

    #[setter(xyzs)]
    pub fn set_xyzs(&mut self, xyzs: Vec<[f64; 3]>) {
        self.core.xyzs = xyzs;
    }

    #[getter(atos)]
    pub fn get_atos(&self) -> Vec<(usize, usize, String, [f64; 3])> {
        self.core.atos.clone()
    }

    #[setter(atos)]
    pub fn set_atos(&mut self, atos: Vec<(usize, usize, String, [f64; 3])>) {
        self.core.atos = atos;
    }

    #[getter(BM)]
    pub fn get_bm(&self) -> Vec<Vec<f64>> {
        // outer_iter 不依赖数组是否连续，也不会像 chunks(ncols) 那样在 ncols==0 时 panic
        self.core.BM.outer_iter().map(|row| row.to_vec()).collect()
    }

    #[setter(BM)]
    pub fn set_bm(&mut self, bm: Vec<Vec<f64>>) -> PyResult<()> {
        self.core.BM = to_array2(&bm, "BM")?;
        Ok(())
    }

    #[getter(CM)]
    pub fn get_cm(&self) -> Vec<Vec<f64>> {
        self.core.CM.outer_iter().map(|row| row.to_vec()).collect()
    }

    #[setter(CM)]
    pub fn set_cm(&mut self, cm: Vec<Vec<f64>>) -> PyResult<()> {
        self.core.CM = to_array2(&cm, "CM")?;
        Ok(())
    }

    #[getter(es)]
    pub fn get_es(&self) -> Vec<f64> {
        self.core.es.clone()
    }

    #[setter(es)]
    pub fn set_es(&mut self, es: Vec<f64>) {
        self.core.es = es;
    }

    pub fn save(&self, path: String) {
        self.core.save(&path);
    }

    pub fn read_mole(&mut self, mole: Mole) {
        self.core.read_mole(&mole.core);
    }
}
