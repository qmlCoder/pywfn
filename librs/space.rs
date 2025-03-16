use rayon::prelude::*;
use std::f64::consts::PI;
use pyo3::prelude::*;
use crate::utils::calc_dist;

// 获得空间格点
fn get_grids(nx:usize,ny:usize,nz:usize) -> Vec<[f64;3]> {
    let mut grids = Vec::with_capacity(nx*ny*nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                grids.push([i as f64,j as f64,k as f64]);
            }
        }
    }
    grids
}

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
    xyzs: &Vec<[f64;3]>, // 原子轨道坐标
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

// 计算一个点处的所有一个分子轨道波函数
fn obt_wfn(
    xyz: &[f64;3], // 该点的坐标
    xyzs: &Vec<[f64;3]>,
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    coefs: &Vec<f64>, // 分子轨道的系数
    level:u32
)->(f64,[f64; 3],[[f64; 3]; 3]){
    let nato = coefs.len(); // 原子轨道数量
    let (ato_wfn0,ato_wfn1,ato_wfn2)=ato_wfn(&xyz, &xyzs, &lmns, &coes, &alps,level);
    let mut obt_wfn0 = 0.0;
    let mut obt_wfn1=[0.0;3];
    let mut obt_wfn2=[[0.0;3];3];
    for ato in 0..nato { // 循环所有的原子轨道
        let coef= coefs[ato];
        if coef==0.0 {continue}; // 如果系数为0，则跳过
        obt_wfn0 += coef * ato_wfn0[ato];
        if level==0 {continue};
        for i in 0..3 {
            obt_wfn1[i] += coef * ato_wfn1[ato][i];
        }
        if level==1 {continue};
        for i in 0..3 {
            for j in 0..3 {
                obt_wfn2[i][j] += coef * ato_wfn2[ato][i][j];
            }
        }
    }
    (obt_wfn0,obt_wfn1,obt_wfn2)
}

// 计算一个点处的电子密度，这些只计算一个点的函数不需要使用numpy数组传参，但是如果有格点的话就需要了
fn mol_rho(
    xyz: &[f64;3],
    xyzs: &Vec<[f64;3]>,
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    matc: &Vec<Vec<f64>>,
    level:u32
) -> (f64,[f64; 3],[[f64; 3]; 3]) {
    let mut mol_rho0 = 0.0;
    let mut mol_rho1=[0.0;3];
    let mut mol_rho2=[[0.0;3];3];
    let nato = matc.len(); //分子轨道数量，传入的系数矩阵要转置了
    let nobt = matc[0].len(); // 原子轨道数量

    let (ato_wfn0,ato_wfn1,ato_wfn2)=ato_wfn(xyz, xyzs, lmns, coes, alps,level);
    let mut obt_wfn0 = 0.0;
    let mut obt_wfn1=[0.0;3];
    let mut obt_wfn2=[[0.0;3];3];
    for obt in 0..nobt {
        for ato in 0..nato { // 循环所有的原子轨道
            let coef= matc[ato][obt];
            if coef == 0.0 {continue}; // 如果系数为0，则跳过
            obt_wfn0 += coef * ato_wfn0[ato];
            if level==0 {continue};
            for i in 0..3 {
                obt_wfn1[i] += coef * ato_wfn1[ato][i];
            }
            if level==1 {continue};
            for i in 0..3 {
                for j in 0..3 {
                    obt_wfn2[i][j] += coef * ato_wfn2[ato][i][j];
                }
            }
        }
        mol_rho0 += obt_wfn0.powf(2.0);
        if level>=1{
            for k in 0..3 {
                mol_rho1[k] += 2.0*obt_wfn0*obt_wfn1[k];
            }
        }
        if level>=2{
            for k in 0..3 {
                for l in 0..3 {
                    mol_rho2[k][l] += 2.0*obt_wfn1[k]*obt_wfn1[l] + 2.0*obt_wfn0*obt_wfn2[k][l];
                }
            }
        }
        obt_wfn0 = 0.0;
        obt_wfn1 = [0.0;3];
        obt_wfn2 = [[0.0;3];3];
    }
    (mol_rho0,mol_rho1,mol_rho2)
}

// 计算所有原子轨道的电子密度
fn ato_rho(
    xyz: &[f64;3],
    xyzs: &Vec<[f64;3]>, // 原子轨道坐标
    lmns: &Vec<[i32;3]>,
    coes: &Vec<Vec<f64>>,
    alps: &Vec<Vec<f64>>,
    matp: &Vec<Vec<f64>>, // 密度矩阵
    level:u32
)->(Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3]; 3]>){
    let nato=xyzs.len();
    let (ato_wfn0,ato_wfn1,ato_wfn2)=ato_wfn(xyz, xyzs, lmns, coes, alps,level);
    let mut ato_rho0=vec![0.0;nato];
    let mut ato_rho1=vec![[0.0;3];nato];
    let mut ato_rho2=vec![[[0.0;3];3];nato];
    for i in 0..nato {
        for j in 0..nato {
            if matp[i][j]==0.0 {continue};
            ato_rho0[i]+=ato_wfn0[i]*ato_wfn0[j]*matp[i][j];
            if level>=1{
                for k in 0..3 {
                    ato_rho1[i][k]+=(ato_wfn0[i]*ato_wfn1[j][k]+ato_wfn0[j]*ato_wfn1[i][k])*matp[i][j];
                }
            }
            if level>=2{
                for k in 0..3 {
                    for l in 0..3 {
                        ato_rho2[i][k][l]+=(
                            ato_wfn2[i][k][l]*ato_wfn0[j]+
                            ato_wfn1[i][k]*ato_wfn1[j][l]+
                            ato_wfn1[i][l]*ato_wfn1[j][k]+
                            ato_wfn0[i]*ato_wfn2[j][k][l])*matp[i][j];
                    }
                }
            }
        }
    }
    (ato_rho0,ato_rho1,ato_rho2)
}

#[pyfunction] // 计算原子中的一个格点再整个分子中的权重
pub fn a2m_weits_rs(
    iatm:usize,
    atm_grids:Vec<[f64;3]>,// 原子格点坐标
    atm_weits:Vec<f64>, // 原子格点权重
    atm_pos:Vec<[f64;3]>, // 原子坐标
    atm_rad:Vec<f64>, // 原子半径
)->PyResult<Vec<f64>>{
    let ngrid=atm_grids.len();
    let natm=atm_pos.len();
    // println!("格点数量：{ngrid}, 原子数量：{natm}");
    // let a2m_grids=vec![[0.0;3];ngrid];
    let mut a2m_weits=Vec::with_capacity(ngrid);
    let mut su: Vec<Vec<f64>>;
    for g in 0..ngrid {
        let gp=atm_grids[g];
        su=vec![vec![1.0;natm];natm];
        for i in 0..natm {
            let pi=atm_pos[i];
            let ri=calc_dist(&pi, &gp);
            for j in 0..natm {
                if j==i {continue};
                let pj=atm_pos[j];
                let rj=calc_dist(&pj, &gp);
                let dij=calc_dist(&pi, &pj);
                let miu_ij=(ri-rj)/dij;
                let chi=atm_rad[i]/atm_rad[j]; // 两原子半径之比
                let mut nu_ij;
                if (chi-1.0).abs()<1e-6 {
                    nu_ij=miu_ij;
                }else{
                    let u_ij=(chi-1.0)/(chi+1.0);
                    let mut a_ij=u_ij/(u_ij.powi(2)-1.0);
                    if a_ij >  0.5 {a_ij=0.5};
                    if a_ij < -0.5 {a_ij=-0.5};
                    nu_ij=miu_ij+a_ij*(1.0-miu_ij.powi(2));
                }
                nu_ij=1.5*nu_ij-0.5*nu_ij.powi(3);
                nu_ij=1.5*nu_ij-0.5*nu_ij.powi(3);
                nu_ij=1.5*nu_ij-0.5*nu_ij.powi(3);

                su[j][i]=0.5*(1.0-nu_ij);
                // println!("{i:>3},{j:>3},{:>10.4}",su[j][i]);
            }
        }
        // 计算权重乘积
        let mut wt = vec![1.0; natm];
        for i in 0..natm {
            for j in 0..natm {
                wt[j] *= su[i][j];
                // println!("是否对称:{}=={}?",su[i][j],su[j][i]);
            }
        }
        // println!("wt:{:>10.4},{:>10.4},{:>10.4}",wt[0],wt[1],wt[2]);
        
        // 归一化并保存结果
        let sum_wt: f64 = wt.iter().sum();
        let rat = if sum_wt.abs() > 1e-12 {
            wt[iatm] / sum_wt
        } else {
            0.0
        };
        
        // a2m_grids.push(gp);
        // println!("{g},swt={sum_wt},rat={rat}");
        a2m_weits.push(atm_weits[g] * rat);
        
    }
    println!("原子到分子权重计算完成");
    Ok(a2m_weits)
}

pub fn nuc_potential(
    qpos: &Vec<[f64;3]>, // query position 请求点（要计算的点）的坐标
    xyzs: &Vec<[f64;3]>,
    nucs: &Vec<f64>,
)->Vec<f64>{ // 计算每个格点的核势能
    let npos=qpos.len();
    let natm=xyzs.len();
    let mut vals=vec![0.0; npos];
    for g in 0..npos {
        for a in 0..natm {
            let dist = calc_dist(&qpos[g], &xyzs[a]); // 计算距离
            if dist<1e-6{continue;}
            vals[g] += nucs[a] / dist;
        }
    }
    vals
}

pub fn ele_potential(
    qpos: &Vec<[f64;3]>, // query position 请求点（要计算的点）的坐标
    grids: &Vec<[f64;3]>,
    weits: &Vec<f64>,
    dens: &Vec<f64>,
)->Vec<f64>{ // 计算每个格点的电子势能
    let npos=qpos.len();
    let ngrid=grids.len();
    let mut vals=vec![0.0; npos];
    for q in 0..npos {
        for g in 0..ngrid {
            let mut dist=calc_dist(&qpos[q], &grids[g]);
            if dist<1e-6 {dist=1e-6;}
            vals[q] += weits[g] * dens[g] / dist;
        }
    }
    vals
}



#[pyfunction]
pub fn obt_wfns_rs(
    grids: Vec<[f64;3]>,
    xyzs: Vec<[f64;3]>,
    lmns: Vec<[i32;3]>,
    coes: Vec<Vec<f64>>,
    alps: Vec<Vec<f64>>,
    coefs: Vec<f64>,
    level:u32
)->PyResult<(Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3];3]>)>{
    let ngrid=grids.len();
    let mut obt_wfn0s=Vec::with_capacity(ngrid);
    let mut obt_wfn1s=Vec::with_capacity(ngrid);
    let mut obt_wfn2s=Vec::with_capacity(ngrid);
    // 并行计算每个网格点的波函数值、梯度和Hessian
    let results: Vec<_> = grids.par_iter()
        .map(|grid| obt_wfn(grid, &xyzs, &lmns, &coes, &alps, &coefs, level))
        .collect();

    // 拆分结果到三个向量
    for (wfn0, wfn1, wfn2) in results {
        obt_wfn0s.push(wfn0);
        obt_wfn1s.push(wfn1);
        obt_wfn2s.push(wfn2);
    }

    Ok((obt_wfn0s, obt_wfn1s, obt_wfn2s))
}

#[pyfunction]
pub fn mol_rhos_rs(
    grids: Vec<[f64;3]>,
    xyzs: Vec<[f64;3]>,
    lmns: Vec<[i32;3]>,
    coes: Vec<Vec<f64>>,
    alps: Vec<Vec<f64>>,
    mat_c: Vec<Vec<f64>>,
    level:u32
)->PyResult<(Vec<f64>,Vec<[f64; 3]>,Vec<[[f64; 3];3]>)>{
    let ngrid=grids.len();
    let mut mol_rho0s=Vec::with_capacity(ngrid);
    let mut mol_rho1s=Vec::with_capacity(ngrid);
    let mut mol_rho2s=Vec::with_capacity(ngrid);
    // 并行计算每个网格点的电子密度、梯度和Hessian
    let results: Vec<_> = grids.par_iter()
        .map(|grid| mol_rho(grid, &xyzs, &lmns, &coes, &alps, &mat_c, level))
        .collect();

    // 将结果拆分到三个向量中
    for (rho0, rho1, rho2) in results {
        mol_rho0s.push(rho0);
        mol_rho1s.push(rho1);
        mol_rho2s.push(rho2);
    }
    Ok((mol_rho0s,mol_rho1s,mol_rho2s))
}

#[pyfunction]
pub fn ato_rhos_rs(
    grids: Vec<[f64;3]>,
    xyzs: Vec<[f64;3]>,
    lmns: Vec<[i32;3]>,
    coes: Vec<Vec<f64>>,
    alps: Vec<Vec<f64>>,
    mat_c: Vec<Vec<f64>>,
    level:u32
)->PyResult<(Vec<Vec<f64>>,Vec<Vec<[f64; 3]>>,Vec<Vec<[[f64; 3]; 3]>>)>{
    let ngrid=grids.len();
    let mut ato_rho0s: Vec<Vec<f64>>=Vec::with_capacity(ngrid);
    let mut ato_rho1s: Vec<Vec<[f64; 3]>>=Vec::with_capacity(ngrid);
    let mut ato_rho2s: Vec<Vec<[[f64; 3]; 3]>>=Vec::with_capacity(ngrid);
    for g in 0..ngrid {
        let (rho0, rho1, rho2) = ato_rho(&grids[g], &xyzs, &lmns, &coes, &alps, &mat_c, level);
        ato_rho0s.push(rho0);
        ato_rho1s.push(rho1);
        ato_rho2s.push(rho2);
    }
    Ok((ato_rho0s,ato_rho1s,ato_rho2s))
}

#[pyfunction]
pub fn nuc_potential_rs(
    qpos: Vec<[f64;3]>,
    xyzs: Vec<[f64;3]>,
    nucs: Vec<f64>,
)->PyResult<Vec<f64>>{
    Ok(nuc_potential(&qpos, &xyzs, &nucs))
}

#[pyfunction]
pub fn ele_potential_rs(
    qpos: Vec<[f64;3]>,
    grids: Vec<[f64;3]>,
    weits: Vec<f64>,
    dens: Vec<f64>,
)->PyResult<Vec<f64>>{
    Ok(ele_potential(&qpos, &grids, &weits, &dens))
}

#[pyfunction]
pub fn get_grids_rs(nx:usize,ny:usize,nz:usize)->PyResult<Vec<[f64;3]>>{
    Ok(get_grids(nx,ny,nz))
}