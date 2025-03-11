pub mod march;
pub mod datas;

use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
use std::f64::consts::PI;

// 计算波函数
// xyz:以原子为中心的格点坐标
fn gtf(
    xyz: &[f64;3],
    lmn: &[i32;3],
    alp: f64,
    level: u32,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let x = xyz[0];
    let y = xyz[1];
    let z = xyz[2];
    let l = lmn[0];
    let m = lmn[1];
    let n = lmn[2];
    let facs = [1., 1., 3.];
    let fac = facs[l as usize] * facs[m as usize] * facs[n as usize];
    let r2 = x.powf(2.0) + y.powf(2.0) + z.powf(2.0);
    let ang = l + m + n;
    let nm = (2.0 * alp / PI).powf(0.75) * ((4.0 * alp).powf(ang as f64) / fac).powf(0.5);

    let wfn0 = nm * x.powf(l as f64) * y.powf(m as f64) * z.powf(n as f64) * (-alp * r2).exp();
    let mut wfn1 = [0.0; 3];
    let mut wfn2 = [[0.0; 3]; 3];
    
    if level >= 1 {
        let dx = if x == 0.0 {
            0.0
        } else {
            (l as f64 / x - 2.0 * alp * x) * wfn0
        };
        let dy = if y == 0.0 {
            0.0
        } else {
            (m as f64 / y - 2.0 * alp * y) * wfn0
        };
        let dz = if z == 0.0 {
            0.0
        } else {
            (n as f64 / z - 2.0 * alp * z) * wfn0
        };
        wfn1 = [dx, dy, dz];
    };
    if level >= 2 {
        let dxx = if x == 0.0 {
            0.0
        } else {
            let t1 = -(l as f64) / x.powi(2) - 2.0 * alp;
            let t2 = (l as f64 / x - 2.0 * alp * x).powi(2);
            (t1 + t2) * wfn0
        };

        let dyy = if y == 0.0 {
            0.0
        } else {
            let t1 = -(m as f64) / y.powi(2) - 2.0 * alp;
            let t2 = (m as f64 / y - 2.0 * alp * y).powi(2);
            (t1 + t2) * wfn0
        };

        let dzz = if z == 0.0 {
            0.0
        } else {
            let t1 = -(n as f64) / z.powi(2) - 2.0 * alp;
            let t2 = (n as f64 / z - 2.0 * alp * z).powi(2);
            (t1 + t2) * wfn0
        };

        let dxy = if x == 0.0 || y == 0.0 {
            0.0
        } else {
            (l as f64 / x - 2.0 * alp * x) * (m as f64 / y - 2.0 * alp * y) * wfn0
        };

        let dxz = if x == 0.0 || z == 0.0 {
            0.0
        } else {
            (l as f64 / x - 2.0 * alp * x) * (n as f64 / z - 2.0 * alp * z) * wfn0
        };

        let dyz = if y == 0.0 || z == 0.0 {
            0.0
        } else {
            (m as f64 / y - 2.0 * alp * y) * (n as f64 / z - 2.0 * alp * z) * wfn0
        };
        wfn2 = [[dxx, dxy, dxz], [dxy, dyy, dyz], [dxz, dyz, dzz]];
    };
    (wfn0, wfn1, wfn2)
}

// 计算收缩基函数（原子轨道），基函数的线性组合
fn cgf(
    xyz: &[f64;3],
    lmn: &[i32;3],
    coes: &Vec<f64>,
    alps: &Vec<f64>,
    level: u32,
) -> (f64, [f64; 3], [[f64; 3]; 3]) {
    let mut wfn0 = 0.0;
    let mut wfn1 = [0.0; 3];
    let mut wfn2 = [[0.0; 3]; 3];
    for i in 0..coes.len() {
        let (val0, val1, val2) = gtf(&xyz, &lmn, alps[i], level);
        wfn0 += coes[i] * val0;
        if level==0 {continue};
        for j in 0..3 {
            wfn1[j] += coes[i] * val1[j];
        }
        if level==1 {continue};
        for j in 0..3 {
            for k in 0..3 {
                wfn2[j][k] += coes[i] * val2[j][k];
            }
        }
    }
    (wfn0, wfn1, wfn2)
}

// 计算一个点处的所有原子轨道的波函数
// xyzs: 原子轨道坐标，coes和alps：每个基函数的系数和指数，每个基函数的数量是不确定的

fn ato_wfn(
    xyz: &[f64;3],
    xyzs: &Vec<[f64;3]>,
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    level:u32
) -> (Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3]; 3]>) {
    let mut wfn0 = Vec::new();
    let mut wfn1 = Vec::new();
    let mut wfn2 = Vec::new();
    for i in 0..xyzs.len() {
        let dx = xyz[0] - xyzs[i][0];
        let dy = xyz[1] - xyzs[i][1];
        let dz = xyz[2] - xyzs[i][2];
        let coes = &coes[i];
        let alps = &alps[i];
        let (val0,val1,val2) = cgf(&[dx, dy, dz], &lmns[i], &coes, &alps,level);
        wfn0.push(val0);
        wfn1.push(val1);
        wfn2.push(val2);
    }
    (wfn0, wfn1, wfn2)
}

// 计算一个点处的所有分子轨道波函数
fn obt_wfn(
    xyz: &[f64;3],
    xyzs: &Vec<[f64;3]>,
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    mat_c: &Vec<Vec<f64>>,
    level:u32
)->(Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3]; 3]>){
    let nato = mat_c.len(); // 原子轨道数量
    let nobt = mat_c[0].len(); //分子轨道数量
    let (ato_wfn0,ato_wfn1,ato_wfn2)=ato_wfn(xyz, xyzs, lmns, coes, alps,level);
    let mut obt_wfn0=Vec::new();
    let mut obt_wfn1=Vec::new();
    let mut obt_wfn2=Vec::new();

    let mut tmp_wfn0 = 0.0;
    let mut tmp_wfn1=[0.0;3];
    let mut tmp_wfn2=[[0.0;3];3];
    for obt in 0..nobt {
        for ato in 0..nato { // 循环所有的原子轨道
            let coef= mat_c[ato][obt];
            tmp_wfn0 += coef * ato_wfn0[ato];
            if level==0 {continue};
            for i in 0..3 {
                tmp_wfn1[i] += coef * ato_wfn1[ato][i];
            }
            if level==1 {continue};
            for i in 0..3 {
                for j in 0..3 {
                    tmp_wfn2[i][j] += coef * ato_wfn2[ato][i][j];
                }
            }
        }
        obt_wfn0.push(tmp_wfn0);
        obt_wfn1.push(tmp_wfn1);
        obt_wfn2.push(tmp_wfn2);
        tmp_wfn0 = 0.0;
        tmp_wfn1=[0.0;3];
        tmp_wfn2=[[0.0;3];3];
    }
    (obt_wfn0, obt_wfn1, obt_wfn2)
}

// 计算一个点处的电子密度，这些只计算一个点的函数不需要使用numpy数组传参，但是如果有格点的话就需要了
fn mol_rho(
    xyz: &[f64;3],
    xyzs: &Vec<[f64;3]>,
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    mat_c: &Vec<Vec<f64>>,
    level:u32
) -> (f64,[f64; 3],[[f64; 3]; 3]) {
    let mut mol_rho0 = 0.0;
    let mut mol_rho1=[0.0;3];
    let mut mol_rho2=[[0.0;3];3];
    let nobt = mat_c[0].len(); //分子轨道数量

    let (obt_wfn0,obt_wfn1,obt_wfn2)=obt_wfn(xyz, xyzs, lmns, coes, alps,mat_c,level);
    for obt in 0..nobt {
        mol_rho0 += obt_wfn0[obt].powf(2.0);
        if level==0 {continue};
        for k in 0..3 {
            mol_rho1[k] += 2.0*obt_wfn0[obt]*obt_wfn1[obt][k];
        }
        if level==1 {continue};
        for k in 0..3 {
            for l in 0..3 {
                mol_rho2[k][l] += 2.0*obt_wfn1[obt][k]*obt_wfn1[obt][l] + 2.0*obt_wfn0[obt]*obt_wfn2[obt][k][l];
            }
        }
    }
    (mol_rho0,mol_rho1,mol_rho2)
}

#[pyfunction]
fn mol_rhos(
    grids: Vec<[f64;3]>,
    xyzs: Vec<[f64;3]>,
    lmns: Vec<[i32;3]>,
    coes: Vec<Vec<f64>>,
    alps: Vec<Vec<f64>>,
    mat_c: Vec<Vec<f64>>,
    level:u32
)->PyResult<(Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3];3]>)>{
    let ngrid=grids.len();
    let mut mol_rho0s=Vec::new();
    let mut mol_rho1s=Vec::new();
    let mut mol_rho2s=Vec::new();
    for i in 0..ngrid {
        let (mol_rho0,mol_rho1,mol_rho2)=mol_rho(&grids[i], &xyzs, &lmns, &coes, &alps,&mat_c,level);
        mol_rho0s.push(mol_rho0);
        mol_rho1s.push(mol_rho1);
        mol_rho2s.push(mol_rho2);
    };
    Ok((mol_rho0s,mol_rho1s,mol_rho2s))
}

#[pymodule]
fn rlib(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(mol_rhos, m)?)?;
    m.add_function(wrap_pyfunction!(march::march_cube, m)?)?;
    Ok(())
}
