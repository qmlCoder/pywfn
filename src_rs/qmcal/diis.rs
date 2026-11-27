use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::prelude::*;
use rswfn;

/// Pulay DIIS（Direct Inversion in the Iterative Subspace）外推器
#[pyclass]
pub struct Diis {
    core: rswfn::qmcal::diis::Diis,
}

#[pymethods]
impl Diis {
    /// 创建 DIIS 实例，`max_len` 为保留的历史向量数。
    #[new]
    pub fn new(max_len: usize) -> Self {
        Self {
            core: rswfn::qmcal::diis::Diis::new(max_len),
        }
    }

    /// 清空历史向量。
    pub fn reset(&mut self) {
        self.core.reset()
    }

    /// 计算 Fock 矩阵与密度矩阵的交换子 `[F, D]`（作为误差矩阵）。
    #[staticmethod]
    pub fn commutator(
        py: Python,
        fock: PyReadonlyArray2<f64>,
        dmat: PyReadonlyArray2<f64>,
        smat: PyReadonlyArray2<f64>,
    ) -> Py<PyArray2<f64>> {
        let fock = fock.as_array().to_owned();
        let dmat = dmat.as_array().to_owned();
        let smat = smat.as_array().to_owned();
        let res = rswfn::qmcal::diis::Diis::commutator(&fock, &dmat, &smat);
        res.into_pyarray(py).unbind()
    }

    /// 用历史向量外推新的 Fock 矩阵；历史不足时返回 `None`。
    pub fn extrapolate(
        &mut self,
        py: Python,
        fock: PyReadonlyArray2<f64>,
        err: PyReadonlyArray2<f64>,
    ) -> Option<Py<PyArray2<f64>>> {
        let fock = fock.as_array().to_owned();
        let err = err.as_array().to_owned();
        self.core
            .extrapolate(&fock, &err)
            .map(|res| res.into_pyarray(py).unbind())
    }
}

/// 在父模块中注册 `diis` 子模块。
pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "diis")?;
    m.add_class::<Diis>()?;
    parent_module.add_submodule(&m)
}
