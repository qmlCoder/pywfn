// 将一些波函数分析的功能直接实现了

use crate::maths::matrix_phalf;
use ndarray::Array2;
// use rayon::vec;

// 计算mulliken电荷，返回电子数量
pub fn mulliken(
    pmat: &Array2<f64>, // 密度矩阵
    smat: &Array2<f64>, // 重叠矩阵
    nums: &Vec<u32>,    // 每个原子基函数数量
) -> Vec<f64> {
    let mut eles = vec![0.0; nums.len()];
    let mut loc = 0;
    let nbas = pmat.nrows(); // 基函数数量
    for a in 0..nums.len() {
        let count = nums[a] as usize;
        for i in loc..loc + count {
            for j in 0..nbas {
                eles[a] += pmat[(i, j)] * smat[(i, j)];
            }
        }
        loc += count;
    }
    eles
}

// 计算lowdin 电荷，返回电子数量
pub fn lowdin(
    pmat: &Array2<f64>, // 密度矩阵
    smat: &Array2<f64>, // 重叠矩阵
    nums: &Vec<u32>,    // 每个原子基函数数量
) -> Vec<f64> {
    let s_half = matrix_phalf(smat); // 重叠矩阵的平方根
    let qmat = pmat.dot(&s_half);

    let natm = nums.len();
    let nmat = smat.shape()[0];
    let mut eles = vec![0.0; nmat];
    let nbas = pmat.nrows();

    let mut loc = 0;
    for a in 0..natm {
        let count = nums[a] as usize;
        for i in loc..loc + count {
            for j in 0..nbas {
                eles[a] += qmat[(i, j)] * s_half[(i, j)];
            }
        }
        loc += count;
    }
    eles
}

// // 计算hirshfeld电荷
// pub fn hirshfeld(
//     atms: &Vec<f64>, // 原子电荷
//     xyzs: &Vec<[f64; 3]>,
// ) {
//     let natm = atms.len();
// }
