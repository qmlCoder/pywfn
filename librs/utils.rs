use pyo3::prelude::*;


pub fn calc_dist(p1:&[f64;3],p2:&[f64;3])->f64 { // 计算两点之间的距离
    let [x0,y0,z0]=p1;
    let [x1,y1,z1]=p2;
    let dist=((x1-x0).powi(2)+(y1-y0).powi(2)+(z1-z0).powi(2)).sqrt();
    dist
}

// 拉格朗日插值
pub fn lag_intpol(xs:&Vec<f64>,ys:&Vec<f64>,ts:&Vec<f64>)->Vec<f64>{
    let nx=xs.len();
    let nt=ts.len();
    let mut vs: Vec<f64>=Vec::with_capacity(nt);
    // let mut loc;
    for t in 0..nt{
        let mut loc=0; // 找到离当前点最近的点
        for i in 0..nx{
            loc=i;
            if xs[i]>ts[t]{break;}
        }
        // println!("{t},loc:{loc}");
        let up=0.max(loc as i32-3) as usize;
        let lo=(nx).min(loc+2);
        // println!("t={t},loc={loc},up={up},lo={lo},{:>10.4}",ts[t]);
        if ts[t]<xs[0]{
            let slop=(ys[1]-ys[0])/(xs[1]-xs[0]);
            vs.push(ys[0]+slop*(ts[t]-xs[0]));
            // println!("0:{}",ys[0]+slop*(ts[t]-xs[0]));
        }else 
        if ts[t]>xs[nx-1] {
            vs.push(0.0);
            // println!("1:{},{}?{}",0.0,ts[t],xs[nx-1]);
        }else {
            let mut tv: f64=0.0;
            for i in up..lo{
                let mut pv=1.0;
                for j in up..lo{
                    if i==j {continue;}
                    pv*=(ts[t]-xs[j])/(xs[i]-xs[j]);
                }
                tv+=pv*ys[i];
                // println!("pv{i:>5}{:>10.4}",tv);
            }
            // println!("2:{}",tv);
            vs.push(tv);
        }
        // println!("{t:>5}{loc:>5}{up:>5}{lo:>5}{:>10.4}{:>10.4}",ts[t],vs[t]);
    }
    vs
}

#[pyfunction]
pub fn lag_intpol_rs(xs:Vec<f64>,ys:Vec<f64>,ts:Vec<f64>)->PyResult<Vec<f64>>{
    let vs=lag_intpol(&xs, &ys, &ts);
    Ok(vs)
}