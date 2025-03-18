use pyo3::PyResult;
use pyo3::prelude::*;

// 此模块定义积分函数
use crate::datas::{rhm,whm};

pub fn gtf_integ(alps:&[f64;2],xyzs:&[[f64;3];2],lmns:&[[i32;3];2])->f64{
    let [e1,e2]=alps;
    let [x1,y1,z1]=&xyzs[0];
    let [x2,y2,z2]=&xyzs[1];
    let alp=e1+e2;
    let px=(e1*x1+e2*x2)/alp;
    let py=(e1*y1+e2*y2)/alp;
    let pz=(e1*z1+e2*z2)/alp;
    let [l1,m1,n1]=lmns[0];
    let [l2,m2,n2]=lmns[1];
    let dx=x2-x1;
    let dy=y2-y1;
    let dz=z2-z1;
    let expv=(-e1*e2*(dx.powi(2)+dy.powi(2)+dz.powi(2))/alp).exp();
    let nx=((l1+l2+1) as f64/2.0).ceil() as usize;
    let ny=((m1+m2+1) as f64/2.0).ceil() as usize;
    let nz=((n1+n2+1) as f64/2.0).ceil() as usize;
    let sqrv=alp.sqrt();
    let rhm=rhm();
    let whm=whm();
    let mut sx=0.0;
    for i in 0..nx{
        let tmpv=rhm[nx-1][i]/sqrv+px;
        sx += whm[nx-1][i]*(tmpv-x1).powi(l1)*(tmpv-x2).powi(l2);
    }
    sx=sx/sqrv;

    let mut sy=0.0;
    for i in 0..ny{
        let tmpv=rhm[ny-1][i]/sqrv+py;
        sy+=whm[ny-1][i]*(tmpv-y1).powi(m1)*(tmpv-y2).powi(m2);
    }
    sy=sy/sqrv;

    let mut sz=0.0;
    for i in 0..nz{
        let tmpv=rhm[nz-1][i]/sqrv+pz;
        sz+=whm[nz-1][i]*(tmpv-z1).powi(n1)*(tmpv-z2).powi(n2);
    }
    sz=sz/sqrv;
    let val=sx*sy*sz*expv;
    // println!("{nx:>3}{ny:>3}{nz:>3}{sx:>10.4}{sy:>10.4}{sz:>10.4}{sqrv:>10.4}{expv:>10.4}");
    val
}

pub fn mat_integ(atos:&Vec<i32>,coes:&Vec<f64>,alps:&Vec<f64>,lmns:&Vec<[i32;3]>,xyzs:&Vec<[f64;3]>)->Vec<Vec<f64>>{
    
    let nato=xyzs.len(); // 原子轨道数量
    let nbas=atos.len(); // 基函数数量
    let mut mats=vec![vec![0.0;nato];nato]; // 重叠矩阵
    let mut bas_ul: Vec<Vec<usize>>=vec![vec![0;2];nato]; //每个原子轨道对应基函数的上下界
    // println!("nato:{nato:>3},nbas:{nbas:>3}");
    let mut iato;
    bas_ul[0][0]=atos[0] as usize;
    for i in 0..nbas{
        // println!("coe,alp{i:>3}{:>10.4}{:>10.4}",coes[i],alps[i]);
        iato=atos[i] as usize;
        if iato==0 {continue;}
        if bas_ul[iato as usize][0]==0{
            bas_ul[iato][0]=i;
            bas_ul[iato-1][1]=i;
        }
    }
    bas_ul[nato-1][1]=nbas;
    let mut xyz_ij=[[0.0;3];2];
    for i in 0..nato{
        for j in 0..nato{
            xyz_ij[0]=xyzs[i];
            xyz_ij[1]=xyzs[j];
            for k in bas_ul[i][0] as usize..bas_ul[i][1] as usize{
                for l in bas_ul[j][0] as usize..bas_ul[j][1] as usize{
                    let sval=gtf_integ(&[alps[k],alps[l]],&xyz_ij,&[lmns[k],lmns[l]]);
                    mats[i][j]+=coes[k]*coes[l]*sval;
                    // println!("{i:>3}{j:>3}{k:>3}{l:>3}{sval:>10.4}{:>10.4}{:>10.4}",coes[k],coes[l]);
                }
            }
        }
    }
    mats
}

#[pyfunction]
pub fn mat_integ_rs(atos:Vec<i32>,coes:Vec<f64>,alps:Vec<f64>,lmns:Vec<[i32;3]>,xyzs:Vec<[f64;3]>)->PyResult<Vec<Vec<f64>>>{
    Ok(mat_integ(&atos,&coes,&alps,&lmns,&xyzs))
}