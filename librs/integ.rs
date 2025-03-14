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
    let nx=((l1+l2+1) as f64/2.0) as usize;
    let ny=((m1+m2+1) as f64/2.0) as usize;
    let nz=((n1+n2+1) as f64/2.0) as usize;
    let sqrv=alp.sqrt();
    let rhm=rhm();
    let whm=whm();
    let mut sx=0.0;
    for i in 0..=nx{
        let tmpv=rhm[nx][i]/sqrv+px;
        sx+=whm[nx][i]*(tmpv-x1).powi(l1)*(tmpv-x2).powi(l2);
    }
    sx=sx/sqrv;

    let mut sy=0.0;
    for i in 0..=ny{
        let tmpv=whm[ny][i]/sqrv+py;
        sy+=whm[ny][i]*(tmpv-y1).powi(m1)*(tmpv-y2).powi(m2);
    }
    sy=sy/sqrv;

    let mut sz=0.0;
    for i in 0..=nz{
        let tmpv=whm[nz][i]/sqrv+pz;
        sz+=whm[nz][i]*(tmpv-z1).powi(n1)*(tmpv-z2).powi(n2);
    }
    sz=sz/sqrv;
    let val=sx*sy*sz*expv;
    val
}

pub fn mat_integ(atos:Vec<i32>,coes:Vec<f64>,alps:Vec<f64>,lmns:Vec<[i32;3]>,xyzs:Vec<[f64;3]>)->Vec<Vec<f64>>{
    let mut mats=vec![vec![0.0;atos.len()];atos.len()];
    let nato=xyzs.len(); // 原子轨道数量
    let nbas=atos.len(); // 基函数数量
    let mut bas_ul=vec![vec![0;2];nbas];
    mats
}