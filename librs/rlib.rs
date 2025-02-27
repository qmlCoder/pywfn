use std::f64::consts::PI;
use pyo3::prelude::*;
use pyo3::wrap_pyfunction;
// use ndarray::{Array,Array2,ArrayView}

// 计算波函数
#[pyfunction]
fn gtf(xyz:Vec<f64>,lmn:Vec<i32>,alp:f64)-> PyResult<f64>{
    let x=xyz[0];
    let y=xyz[1];
    let z=xyz[2];
    let l=lmn[0];
    let m=lmn[1];
    let n=lmn[2];
    let face = [1.,2.,3.];
    let fac = face[l as usize]*face[m as usize]*face[n as usize];
    let r2 = x.powf(2.0)+y.powf(2.0)+z.powf(2.0);
    let ang = l+m+n;
    let nm = (2.0*alp/PI).powf(0.75)*((4.0*alp).powf(ang as f64)/fac).powf(0.5);
    
    let wfn0 = nm*x.powf(l as f64)*y.powf(m as f64)*z.powf(n as f64)*(-alp*r2).exp();
    Ok(wfn0)
}

// 计算收缩基函数（原子轨道），基函数的线性组合
#[pyfunction]
fn cgf(xyz:Vec<f64>,lmn:Vec<i32>,coes:Vec<f64>,alps:Vec<f64>)->PyResult<f64>{
    let mut wfn = 0.0;
    for i in 0..coes.len(){
        wfn += coes[i]*gtf(xyz.clone(),lmn.clone(),alps[i])?;
    }
    Ok(wfn)
}

// 计算一个点处的所有原子轨道的波函数
// xyzs: 原子轨道坐标，coes和alps：每个基函数的系数和指数，每个基函数的数量是不确定的
#[pyfunction]
fn ato_wfn(xyz:Vec<f64>,xyzs:Vec<Vec<f64>>,lmns:Vec<Vec<i32>>,coes:Vec<Vec<f64>>,alps:Vec<Vec<f64>>)->PyResult<Vec<f64>>{
    let mut wfns = Vec::new();
    for i in 0..xyzs.len(){
        let dx= xyz[0]-xyzs[i][0];
        let dy= xyz[1]-xyzs[i][1];
        let dz= xyz[2]-xyzs[i][2];
        let coes=coes[i].clone();
        let alps=alps[i].clone();
        wfns.push(cgf([dx,dy,dz].to_vec(),lmns[i].clone(),coes,alps)?);
    }
    Ok(wfns)
}

// 计算一个点处的电子密度，这些只计算一个点的函数不需要使用numpy数组传参，但是如果有格点的话就需要了
#[pyfunction]
fn mol_rho(xyz:Vec<f64>,xyzs:Vec<Vec<f64>>,lmns:Vec<Vec<i32>>,coes:Vec<Vec<f64>>,alps:Vec<Vec<f64>>,mat_c:Vec<Vec<f64>>)->PyResult<f64>{
    let mut mol_rho = 0.0;
    let ato_wfns = ato_wfn(xyz,xyzs,lmns,coes,alps)?;
    let nato = mat_c.len(); // 原子轨道数量
    let nobt = mat_c[0].len(); //分子轨道数量
    println!("nato:{},nobt:{}",nato,nobt);
    let mut obt_rho=0.0;
    for i in 0..nobt{
        for j in 0..nato{   
            obt_rho += mat_c[j][i]*ato_wfns[j];
            // println!("{},{},{},{}",i,j,nato,nobt)
        }
        mol_rho += obt_rho.powf(2.0);
    }
    Ok(mol_rho)
}

#[pymodule]
fn rlib(m: &Bound<'_,PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(gtf,m)?)?;
    m.add_function(wrap_pyfunction!(cgf,m)?)?;
    m.add_function(wrap_pyfunction!(ato_wfn,m)?)?;
    m.add_function(wrap_pyfunction!(mol_rho,m)?)?;
    Ok(())
}