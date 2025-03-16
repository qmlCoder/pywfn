use pyo3::prelude::*;

pub fn ele_mat(matc:&Vec<Vec<f64>>,mats:&Vec<Vec<f64>>)->Vec<Vec<f64>>{
    let nato=matc.len();
    let nobt=matc[0].len();
    let mut matn=vec![vec![0.0; nobt]; nato];
    for a in 0..nato{
        for o in 0..nobt{
            for i in 0..nato{
                matn[a][o] += matc[a][o]*matc[i][o]*mats[i][a];
            }
        }
    }
    println!("nato:{nato},nobt:{nobt}");
    matn
}


#[pyfunction]
pub fn ele_mat_rs(matc:Vec<Vec<f64>>,mats:Vec<Vec<f64>>)->Vec<Vec<f64>>{
    ele_mat(&matc,&mats)
}
