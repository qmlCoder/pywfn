use crate::base::Mole;
use ndarray::{Array1, Array2, Array3};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArray3, PyArrayMethods};
use pyo3::prelude::*;
use rswfn::{self, props::gridprop::GridVals};

fn grids2array3(py: Python, grids: Vec<[f64; 3]>) -> Py<PyArray2<f64>> {
    let ngrid = grids.len();
    Array2::from_shape_vec((ngrid, 3), grids.into_flattened())
        .unwrap()
        .into_pyarray(py)
        .unbind()
}

#[pyclass]
pub struct Calculator {
    pub mole: Mole,
}

impl Calculator {
    /// 创建底层的 `rswfn::props::gridprop::dftgrid::Calculator` 实例。
    fn caler(&self) -> rswfn::props::gridprop::dftgrid::Calculator<'_> {
        rswfn::props::gridprop::dftgrid::Calculator::new(&self.mole.core)
    }
}

#[pymethods]
impl Calculator {
    #[new]
    /// 创建一个新的 DFT 格点计算器。
    ///
    /// * `mole` - 分子对象
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    pub fn atm_grids(
        &self,
        py: Python,
        iatm: usize,
        nrad: usize,
        fsph: usize,
    ) -> PyResult<(Py<PyArray2<f64>>, Py<PyArray1<f64>>)> {
        let (grids, weits) = self.caler().atm_grids(iatm, nrad, fsph);
        Ok((grids2array3(py, grids), weits.into_pyarray(py).unbind()))
    }

    pub fn a2m_weits(
        &self,
        py: Python,
        iatm: usize,
        atm_grids: Vec<[f64; 3]>,
        atm_weits: Vec<f64>,
    ) -> Py<PyArray1<f64>> {
        self.caler()
            .a2m_weits(iatm, &atm_grids, &atm_weits)
            .into_pyarray(py)
            .unbind()
    }

    pub fn mol_grids(
        &self,
        py: Python,
        nrad: usize,
        fsph: usize,
    ) -> PyResult<(Py<PyArray2<f64>>, Py<PyArray1<f64>>)> {
        let (grids, weits) = self.caler().mol_grids(nrad, fsph);
        Ok((grids2array3(py, grids), weits.into_pyarray(py).unbind()))
    }

    /// 计算分子片段（指定原子集合）的 DFT 格点及权重。
    ///
    /// * `frag` - 片段包含的原子索引列表
    ///
    /// 返回 `(片段格点坐标, 片段格点权重)` 二元组。
    pub fn frag_grids(
        &self,
        py: Python,
        frag: Vec<usize>,
    ) -> PyResult<(Py<PyArray2<f64>>, Py<PyArray1<f64>>)> {
        let (grids, weits) = self.caler().frag_grids(&frag);
        Ok((grids2array3(py, grids), weits.into_pyarray(py).unbind()))
    }
}

#[pyfunction]
pub fn rad_grids(
    py: Python,
    R: f64,
    nrad: usize,
) -> PyResult<(Py<PyArray1<f64>>, Py<PyArray1<f64>>)> {
    let (grids, weits) = rswfn::props::gridprop::dftgrid::rad_grids(R, nrad);
    Ok((
        grids.into_pyarray(py).unbind(),
        weits.into_pyarray(py).unbind(),
    ))
}

#[pyfunction]
pub fn sph_grids(py: Python, nsph: usize) -> PyResult<(Py<PyArray2<f64>>, Py<PyArray1<f64>>)> {
    let (grids, weits) = rswfn::props::gridprop::dftgrid::sph_grids(nsph);
    let ngrid = grids.len();
    Ok((grids2array3(py, grids), weits.into_pyarray(py).unbind()))
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "dftgrid")?;
    m.add_class::<Calculator>()?;
    m.add_function(wrap_pyfunction!(rad_grids, &m)?)?;
    m.add_function(wrap_pyfunction!(sph_grids, &m)?)?;
    parent_module.add_submodule(&m)
}
