use std::collections::HashMap;

use crate::reads::log;
use ndarray::Array2;
use pyo3::prelude::*;

fn array2_vector(array: &Array2<f64>) -> Vec<Vec<f64>> {
    let shape = array.shape();
    let nrow = shape[0];
    let ncol = shape[1];
    let mut vec: Vec<Vec<f64>> = Vec::with_capacity(nrow);
    for i in 0..nrow {
        let mut row_vec = Vec::with_capacity(ncol);
        for j in 0..ncol {
            row_vec.push(array[[i, j]]);
        }
        vec.push(row_vec);
    }
    vec
}

#[pyfunction]
pub fn search_title(path: String) -> PyResult<HashMap<String, u32>> {
    let titles = log::search_title(&path);
    Ok(titles)
}

#[pyfunction]
pub fn read_geome(path: String, title_number: u32) -> PyResult<(Vec<u32>, Vec<[f64; 3]>)> {
    let (atms, xyzs) = log::read_geome(&path, title_number);
    Ok((atms, xyzs))
}

#[pyfunction]
pub fn read_nbas(path: String, title_number: u32) -> PyResult<u32> {
    let nbas = log::read_nbas(&path, title_number);
    // println!("{:?}",nbas);
    Ok(nbas)
}

#[pyfunction]
pub fn read_cmat(
    path: String,
    title_number: u32,
    nbas: u32,
) -> PyResult<(
    Vec<u32>,
    Vec<u32>,
    Vec<String>,
    Vec<f64>,
    Vec<bool>,
    Vec<Vec<f64>>,
)> {
    let (ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat) =
        log::read_cmat(&path, title_number, nbas);
    let cmat = array2_vector(&cmat);
    Ok((ato_atms, ato_shls, ato_syms, obt_engs, obt_occs, cmat))
}

#[pyfunction]
pub fn read_basis(path: String, title_number: u32) -> PyResult<Vec<(u32, i32, u32, f64, f64)>> {
    let basis = log::read_basis(&path, title_number);
    Ok(basis)
}
