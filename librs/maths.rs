use ndarray::{Array, Array2};
use ndarray_linalg::{Eigh, UPLO};

/// 计算(重叠)矩阵的平方根 S^{1/2}
pub fn matrix_phalf(smat: &Array2<f64>) -> Array2<f64> {
    let (eig_val, vec_mat) = smat.eigh(UPLO::Lower).expect("矩阵对角化失败");
    // let inv_mat= vec_mat.inv().expect("矩阵求逆失败");
    let nmat = smat.shape()[0];
    let mut val_mat = Array::<f64, _>::zeros((nmat, nmat));

    for i in 0..nmat {
        val_mat[(i, i)] = eig_val[i].sqrt();
    }
    vec_mat.dot(&val_mat).dot(&vec_mat.t())
}

pub fn matrix_nhalf(smat: &Array2<f64>) -> Array2<f64> {
    let (eig_val, vec_mat) = smat.eigh(UPLO::Lower).expect("矩阵对角化失败");
    // let inv_mat= vec_mat.inv().expect("矩阵求逆失败");
    let nmat = smat.shape()[0];
    let mut val_mat = Array::<f64, _>::zeros((nmat, nmat));
    for i in 0..nmat {
        val_mat[(i, i)] = 1.0 / eig_val[i].sqrt();
    }
    vec_mat.dot(&val_mat).dot(&vec_mat.t())
}
