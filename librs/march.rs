
use crate::datas::MARCH;
use pyo3::prelude::*;


#[pyfunction]
pub fn march_cube(shape: [usize; 3], grids: Vec<[f64;3]>, values: Vec<f64>,isov:f64)->PyResult<Vec<[f64;3]>>{
    println!("Marching Cube算法");
    let voxel=grids2voxel(shape, &grids, &values);
    let verts=voxel2verts(&voxel, isov);
    Ok(verts)
}


fn grids2voxel(shape: [usize; 3], grids: &Vec<[f64;3]>, values: &Vec<f64>) -> Vec<Vec<Vec<[f64;4]>>> {
    let [nx,ny,nz]=shape;
    let mut voxel:Vec<Vec<Vec<[f64;4]>>> = vec![vec![vec![[0.0,0.0,0.0,0.0];nx];ny];nz];
    let mut n=0;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let [x,y,z]=grids[n];
                let v=values[n];
                voxel[i][j][k]=[x,y,z,v];
                n+=1;
            }
        }
    }
    println!("格点转体素数据");
    voxel
}


fn voxel2verts(voxel:&Vec<Vec<Vec<[f64;4]>>>,isov:f64)-> Vec<[f64;3]> {
    let nx=voxel.len();
    let ny=voxel[0].len();
    let nz=voxel[0][0].len();
    let mut verts:Vec<[f64;3]> = Vec::new();
    for i in 0..nx-1 {
        for j in 0..ny-1 {
            for k in 0..nz-1 {
                let (xyzs,vals)=get_around(&voxel,i,j,k);
                // 获取vals的最小值及最小值所在的索引
                let vmin=vals.iter().min_by(|a,b|a.partial_cmp(b).unwrap()).unwrap();
                let vmax=vals.iter().max_by(|a,b|a.partial_cmp(b).unwrap()).unwrap();

                let mut key:usize=0;
                // println!("vmin:{vmin} vmax:{vmax} isov:{isov}");
                if *vmin<isov && *vmax>isov {
                    for n in 0..8 {
                        key+=if vals[n]>isov {(2 as usize).pow((7-n) as u32)} else {0};
                        
                    }
                    println!("key:{key},{vals:?}");
                    let box_verts=get_boxverts(key, xyzs, vals, isov);
                    verts.extend(box_verts);
                }
            }
        }
    }
    verts
}

fn get_around(voxel:&Vec<Vec<Vec<[f64;4]>>>,i:usize,j:usize,k:usize)->([[f64;3];8],[f64;8]) {
    let mut xyzs=[[0.0;3];8];
    let mut vals=[0.0;8];
    let didx:[[usize;3];8]=[
        [0,0,0],
        [1,0,0],
        [1,1,0],
        [0,1,0],
        [0,0,1],
        [1,0,1],
        [1,1,1],
        [0,1,1],
    ];
    for n in 0..8 {
        let [di,dj,dk]=didx[n];
        let x=voxel[i+di][j+dj][k+dk][0];
        let y=voxel[i+di][j+dj][k+dk][1];
        let z=voxel[i+di][j+dj][k+dk][2];
        let v=voxel[i+di][j+dj][k+dk][3];
        xyzs[n]=[x,y,z];
        vals[n]=v;
        // println!("{n}:x={x},y={y},z={z},v={v}")
    }
    (xyzs,vals)
}

fn get_boxverts(key:usize,xyzs:[[f64;3];8],vals:[f64;8],isov:f64)->Vec<[f64;3]> {
    let bonds=[[0, 1], [1, 2], [2, 3], [0, 3], [1, 5], [2, 6], [3, 7], [0, 4], [4, 5], [5, 6], [6, 7], [7, 4]];
    let mut verts:Vec<[f64;3]> = Vec::new();
    let nface=MARCH[key].len()/3; // 面的数量，每个面3个顶点
    println!("key:{key} nface:{nface}");
    for f in 0..nface {
        let face=MARCH[key][f*3..(f+1)*3].to_vec();
        println!("face:{f} {face:?}");
        for p in 0..3 {
            let e=face[p];
            let a=bonds[e as usize][0];
            let b=bonds[e as usize][1];
            let pa=xyzs[a as usize];
            let pb=xyzs[b as usize];
            let va=vals[a as usize];
            let vb=vals[b as usize];
            let t=(isov-va)/(vb-va);
            let x=pa[0]+t*(pb[0]-pa[0]);
            let y=pa[1]+t*(pb[1]-pa[1]);
            let z=pa[2]+t*(pb[2]-pa[2]);
            verts.push([x,y,z]);
        }
    }
    verts
}