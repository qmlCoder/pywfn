#![allow(non_snake_case)]

use numpy::{IntoPyArray, PyArray4, PyReadonlyArray4};
use pyo3::prelude::*;
use rswfn;

use crate::base::Vec3;

/// 把 Python 回调包装成 Rust 闭包：输入坐标，返回函数值
fn wrap_func(func: Py<PyAny>) -> impl Fn(&[f64; 3]) -> f64 + Send + Sync {
    move |xyz: &[f64; 3]| {
        Python::attach(|py| {
            func.call1(py, (xyz.to_vec(),))
                .and_then(|res| res.extract::<f64>(py))
                .expect("Python 回调执行失败")
        })
    }
}

/// 立方体区域，用于等值面提取
#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct Cube {
    pub p0: [f64; 3],
    pub p1: [f64; 3],
}

impl Cube {
    /// 按对角点重建 rswfn 的立方体对象
    pub fn core(&self) -> rswfn::maths::march::Cube {
        rswfn::maths::march::Cube::new(self.p0, self.p1)
    }

    pub fn wrap(core: rswfn::maths::march::Cube) -> Self {
        let p0 = core.p0();
        let p1 = core.p1();
        Self {
            p0: p0.into_arr(),
            p1: p1.into_arr(),
        }
    }
}

#[pymethods]
impl Cube {
    /// 通过体对角线上的两个对角点创建立方体
    ///
    /// ## 参数
    /// - `p0`: 起点坐标
    /// - `p1`: 终点坐标
    #[new]
    pub fn new(p0: [f64; 3], p1: [f64; 3]) -> Self {
        Self { p0, p1 }
    }

    /// 按格点数量计算步长
    ///
    /// ## 参数
    /// - `count`: 每个方向的格点数量
    pub fn calc_step(&self, count: usize) -> [f64; 3] {
        self.core().calc_step(count)
    }

    /// 立方体体积
    pub fn volum(&self) -> f64 {
        self.core().volum()
    }

    /// 起点坐标
    pub fn p0(&self) -> Vec3 {
        Vec3 {
            core: self.core().p0(),
        }
    }

    /// 终点坐标
    pub fn p1(&self) -> Vec3 {
        Vec3 {
            core: self.core().p1(),
        }
    }

    /// 按格点数量生成格点
    ///
    /// ## 参数
    /// - `count`: 每个方向的格点数量
    ///
    /// ## 返回
    /// 返回 (形状, 格点坐标列表, 步长)
    pub fn grids_v1(&self, count: usize) -> ([usize; 3], Vec<[f64; 3]>, f64) {
        self.core().grids_v1(count)
    }

    /// 按步长生成格点
    ///
    /// ## 参数
    /// - `step`: 步长
    ///
    /// ## 返回
    /// 返回 (形状, 格点坐标列表)
    pub fn grids_v2(&self, step: f64) -> ([usize; 3], Vec<[f64; 3]>) {
        self.core().grids_v2(step)
    }
}

/// 从三维标量场中提取等值面区域
///
/// ## 参数
/// - `func`: 计算函数值的 Python 回调，接收坐标列表，返回浮点数
/// - `cube`: 立方体区域
/// - `isov`: 等值面阈值
/// - `step`: 步长
/// - `xyzs`: 原子坐标列表，可选
/// - `limit`: 距离限制，可选
/// - `extract`: 提取函数，接收 (形状, 格点, 数值, 阈值) 并返回 (顶点, 类型)
///
/// ## 返回
/// 返回 (顶点坐标, 顶点类型)
#[pyfunction]
#[pyo3(signature = (func, cube, isov, step, xyzs=None, limit=None, extract=None))]
pub fn extract_region(
    py: Python,
    func: Py<PyAny>,
    cube: Cube,
    isov: f64,
    step: f64,
    xyzs: Option<Vec<[f64; 3]>>,
    limit: Option<f64>,
    extract: Option<Py<PyAny>>,
) -> (Vec<[f64; 3]>, Vec<isize>) {
    let f = wrap_func(func);
    let ex = extract;
    let res = py.detach(move || {
        let extract_fn = |shape: &[usize; 3], grids: &Vec<[f64; 3]>, values: &Vec<f64>, iso: f64| {
            match &ex {
                Some(call) => Python::attach(|py| {
                    call.call1(
                        py,
                        (
                            shape.to_vec(),
                            grids.clone(),
                            values.clone(),
                            iso,
                        ),
                    )
                    .and_then(|res| res.extract::<(Vec<[f64; 3]>, Vec<isize>)>(py))
                    .expect("Python 提取回调执行失败")
                }),
                None => (vec![], vec![]),
            }
        };
        rswfn::maths::march::extract_region(
            &f,
            cube.core(),
            isov,
            step,
            xyzs.as_ref(),
            limit,
            extract_fn,
        )
    });
    res
}

/// 把顶点和三角面片展平为一维数组
///
/// ## 参数
/// - `verts`: 顶点坐标列表
/// - `faces`: 三角面片索引列表
///
/// ## 返回
/// 返回 (展平后的顶点, 展平后的面片)
#[pyfunction]
pub fn flatten(verts: Vec<[f64; 3]>, faces: Vec<[usize; 3]>) -> (Vec<f64>, Vec<usize>) {
    rswfn::maths::march::flatten(&verts, &faces)
}

/// 把格点和数值整理成四维体素数组
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `shape`: 格点形状 [nx, ny, nz]
/// - `grids`: 格点坐标列表
/// - `values`: 格点数值列表
///
/// ## 返回
/// 返回形状为 (nx, ny, nz, 4) 的数组
#[pyfunction]
pub fn grids2voxel(
    py: Python,
    shape: [usize; 3],
    grids: Vec<[f64; 3]>,
    values: Vec<f64>,
) -> Py<PyArray4<f64>> {
    rswfn::maths::march::grids2voxel(&shape, &grids, &values)
        .into_pyarray(py)
        .unbind()
}

/// 合并重复顶点
///
/// ## 参数
/// - `old_verts`: 原始顶点坐标列表
///
/// ## 返回
/// 返回 (去重后的顶点, 三角面片)
#[pyfunction]
pub fn merge_verts(old_verts: Vec<[f64; 3]>) -> (Vec<[f64; 3]>, Vec<[usize; 3]>) {
    rswfn::maths::march::merge_verts(&old_verts)
}

/// 平滑等值面顶点
///
/// ## 参数
/// - `verts`: 顶点坐标列表
/// - `types`: 顶点类型列表
/// - `isov`: 等值面阈值
/// - `func`: 计算函数值的 Python 回调
/// - `iterations`: 迭代次数
///
/// ## 返回
/// 返回平滑后的顶点坐标
#[pyfunction]
pub fn smooth_isosurface(
    py: Python,
    verts: Vec<[f64; 3]>,
    types: Vec<isize>,
    isov: f64,
    func: Py<PyAny>,
    iterations: usize,
) -> Vec<[f64; 3]> {
    let f = wrap_func(func);
    py.detach(move || rswfn::maths::march::smooth_isosurface(&verts, &types, isov, &f, iterations))
}

/// 把顶点吸附到等值面上
///
/// ## 参数
/// - `verts`: 顶点坐标列表
/// - `types`: 顶点类型列表
/// - `isov`: 等值面阈值
/// - `func`: 计算函数值的 Python 回调
/// - `h`: 步长
///
/// ## 返回
/// 返回吸附后的顶点坐标
#[pyfunction]
pub fn snap_to_isosurface(
    py: Python,
    mut verts: Vec<[f64; 3]>,
    types: Vec<isize>,
    isov: f64,
    func: Py<PyAny>,
    h: f64,
) -> Vec<[f64; 3]> {
    let f = wrap_func(func);
    py.detach(move || {
        rswfn::maths::march::snap_to_isosurface(&mut verts, &types, isov, &f, h);
        verts
    })
}

/// 按数值合并两组网格
///
/// ## 参数
/// - `pverts`: 第一组顶点
/// - `pfaces`: 第一组面片
/// - `nverts`: 第二组顶点
/// - `nfaces`: 第二组面片
///
/// ## 返回
/// 返回 (顶点, 面片, 类型)
#[pyfunction]
pub fn value_merge(
    pverts: Vec<[f64; 3]>,
    pfaces: Vec<[usize; 3]>,
    nverts: Vec<[f64; 3]>,
    nfaces: Vec<[usize; 3]>,
) -> (Vec<[f64; 3]>, Vec<[usize; 3]>, Vec<isize>) {
    rswfn::maths::march::value_merge(&pverts, &pfaces, &nverts, &nfaces)
}

/// 从体素数组中提取所有包含等值面的立方体
///
/// ## 参数
/// - `voxel`: 四维体素数组
/// - `isov`: 等值面阈值
///
/// ## 返回
/// 返回立方体列表
#[pyfunction]
pub fn voxel2cubes(voxel: PyReadonlyArray4<f64>, isov: f64) -> Vec<Cube> {
    let voxel = voxel.as_array().to_owned();
    rswfn::maths::march::voxel2cubes(&voxel, isov)
        .into_iter()
        .map(Cube::wrap)
        .collect()
}

pub mod mc {
    use super::{wrap_func, Cube};
    use numpy::PyReadonlyArray4;
    use pyo3::prelude::*;
    use rswfn;

    /// Marching Cubes 的三角形连接表
    #[pyfunction]
    fn MARCH() -> Vec<Vec<i8>> {
        rswfn::maths::march::mc::MARCH.iter().map(|row| row.to_vec()).collect()
    }

    /// 使用 Marching Cubes 算法提取等值面
    ///
    /// ## 参数
    /// - `shape`: 三个方向上的格点数量 [nx, ny, nz]
    /// - `grids`: 所有格点的坐标列表
    /// - `values`: 所有格点上的标量值列表
    /// - `isov`: 等值面的值
    ///
    /// ## 返回
    /// 返回 (顶点坐标, 顶点类型)
    #[pyfunction]
    pub fn march_cube_v1(
        shape: [usize; 3],
        grids: Vec<[f64; 3]>,
        values: Vec<f64>,
        isov: f64,
    ) -> (Vec<[f64; 3]>, Vec<isize>) {
        rswfn::maths::march::mc::march_cube_v1(&shape, &grids, &values, isov)
    }

    /// 使用 Marching Cubes 算法提取等值面（直接给定函数）
    ///
    /// ## 参数
    /// - `func`: 计算函数值的 Python 回调
    /// - `cube`: 立方体区域
    /// - `isov`: 等值面阈值
    /// - `step`: 步长
    /// - `xyzs`: 原子坐标列表，可选
    /// - `limit`: 距离限制，可选
    ///
    /// ## 返回
    /// 返回 (顶点坐标, 顶点类型)
    #[pyfunction]
    #[pyo3(signature = (func, cube, isov, step, xyzs=None, limit=None))]
    pub fn march_cube_v2(
        py: Python,
        func: Py<PyAny>,
        cube: Cube,
        isov: f64,
        step: f64,
        xyzs: Option<Vec<[f64; 3]>>,
        limit: Option<f64>,
    ) -> (Vec<[f64; 3]>, Vec<isize>) {
        let f = wrap_func(func);
        py.detach(move || {
            rswfn::maths::march::mc::march_cube_v2(f, cube.core(), isov, step, xyzs.as_ref(), limit)
        })
    }

    /// 提取体素中所有格点及其数值
    ///
    /// ## 参数
    /// - `voxel`: 四维体素数组
    /// - `isov`: 等值面阈值
    ///
    /// ## 返回
    /// 返回 (所有顶点, 活跃顶点)
    #[pyfunction]
    pub fn voxel2verts(
        voxel: PyReadonlyArray4<f64>,
        isov: f64,
    ) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
        let voxel = voxel.as_array().to_owned();
        rswfn::maths::march::mc::voxel2verts(&voxel, isov)
    }

    pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
        let m = PyModule::new(parent_module.py(), "mc")?;
        m.add_function(wrap_pyfunction!(MARCH, &m)?)?;
        m.add_function(wrap_pyfunction!(march_cube_v1, &m)?)?;
        m.add_function(wrap_pyfunction!(march_cube_v2, &m)?)?;
        m.add_function(wrap_pyfunction!(voxel2verts, &m)?)?;
        parent_module.add_submodule(&m)
    }
}

pub mod mt {
    use super::{wrap_func, Cube};
    use pyo3::prelude::*;
    use rswfn;

    /// 使用 Marching Tetrahedra 算法提取等值面
    ///
    /// ## 参数
    /// - `shape`: 三个方向上的格点数量 [nx, ny, nz]
    /// - `grids`: 所有格点的坐标列表
    /// - `values`: 所有格点上的标量值列表
    /// - `isov`: 等值面的值
    ///
    /// ## 返回
    /// 返回 (顶点坐标, 顶点类型)
    #[pyfunction]
    pub fn march_tetra_v1(
        shape: [usize; 3],
        grids: Vec<[f64; 3]>,
        values: Vec<f64>,
        isov: f64,
    ) -> (Vec<[f64; 3]>, Vec<isize>) {
        rswfn::maths::march::mt::march_tetra_v1(&shape, &grids, &values, isov)
    }

    /// 使用 Marching Tetrahedra 算法提取等值面（直接给定函数）
    ///
    /// ## 参数
    /// - `func`: 计算函数值的 Python 回调
    /// - `cube`: 立方体区域
    /// - `isov`: 等值面阈值
    /// - `step`: 步长
    /// - `xyzs`: 原子坐标列表，可选
    /// - `limit`: 距离限制，可选
    ///
    /// ## 返回
    /// 返回 (顶点坐标, 顶点类型)
    #[pyfunction]
    #[pyo3(signature = (func, cube, isov, step, xyzs=None, limit=None))]
    pub fn march_tetra_v2(
        py: Python,
        func: Py<PyAny>,
        cube: Cube,
        isov: f64,
        step: f64,
        xyzs: Option<Vec<[f64; 3]>>,
        limit: Option<f64>,
    ) -> (Vec<[f64; 3]>, Vec<isize>) {
        let f = wrap_func(func);
        py.detach(move || {
            rswfn::maths::march::mt::march_tetra_v2(f, cube.core(), isov, step, xyzs.as_ref(), limit)
        })
    }

    pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
        let m = PyModule::new(parent_module.py(), "mt")?;
        m.add_function(wrap_pyfunction!(march_tetra_v1, &m)?)?;
        m.add_function(wrap_pyfunction!(march_tetra_v2, &m)?)?;
        parent_module.add_submodule(&m)
    }
}

pub mod sn {
    use pyo3::prelude::*;
    use rswfn;

    /// 使用 Surface Nets 算法提取等值面
    ///
    /// ## 参数
    /// - `shape`: 三个方向上的格点数量 [nx, ny, nz]
    /// - `grids`: 所有格点的坐标列表
    /// - `values`: 所有格点上的标量值列表
    /// - `isov`: 等值面的值
    ///
    /// ## 返回
    /// 返回 (顶点坐标, 顶点类型)
    #[pyfunction]
    pub fn surface_nets_v1(
        shape: [usize; 3],
        grids: Vec<[f64; 3]>,
        values: Vec<f64>,
        isov: f64,
    ) -> (Vec<[f64; 3]>, Vec<isize>) {
        rswfn::maths::march::sn::surface_nets_v1(&shape, &grids, &values, isov)
    }

    pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
        let m = PyModule::new(parent_module.py(), "sn")?;
        m.add_function(wrap_pyfunction!(surface_nets_v1, &m)?)?;
        parent_module.add_submodule(&m)
    }
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "march")?; // 为当前的rs文件创建一个子模块
    m.add_class::<Cube>()?;
    m.add_function(wrap_pyfunction!(extract_region, &m)?)?;
    m.add_function(wrap_pyfunction!(flatten, &m)?)?;
    m.add_function(wrap_pyfunction!(grids2voxel, &m)?)?;
    m.add_function(wrap_pyfunction!(merge_verts, &m)?)?;
    m.add_function(wrap_pyfunction!(smooth_isosurface, &m)?)?;
    m.add_function(wrap_pyfunction!(snap_to_isosurface, &m)?)?;
    m.add_function(wrap_pyfunction!(value_merge, &m)?)?;
    m.add_function(wrap_pyfunction!(voxel2cubes, &m)?)?;
    m.add_function(wrap_pyfunction!(march_cube_v1, &m)?)?;
    mc::register_module(&m)?;
    mt::register_module(&m)?;
    sn::register_module(&m)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}

/// 保持旧接口可用：maths.march.march_cube_v1
#[pyfunction]
pub fn march_cube_v1(
    shape: [usize; 3],
    grids: Vec<[f64; 3]>,
    values: Vec<f64>,
    isov: f64,
) -> (Vec<[f64; 3]>, Vec<isize>) {
    rswfn::maths::march::mc::march_cube_v1(&shape, &grids, &values, isov)
}

