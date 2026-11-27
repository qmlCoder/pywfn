#![allow(non_snake_case)]

use std::collections::HashMap;

use crate::base::{Mole, Stm};
use ndarray::Array2;
use numpy::{IntoPyArray, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::prelude::*;
use rswfn;

#[pyclass]
pub struct RsInt {
    pub mole: Mole,
}

impl RsInt {
    // 提取公共的 calculator 创建逻辑
    fn caler(&self) -> rswfn::rsint::RsInt<'_> {
        rswfn::rsint::RsInt::new(&self.mole.core)
    }
}

#[pymethods]
impl RsInt {
    /// 创建一个新的 RsInt 积分计算实例
    ///
    /// ## 参数
    /// - `mole`: 包含分子、基组信息的 Mole 对象
    ///
    /// ## 返回
    /// 返回一个 RsInt 对象
    #[new]
    pub fn new(mole: Mole) -> Self {
        Self { mole }
    }

    /// 计算重叠矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph"（球谐）或 "car"（笛卡尔）
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的重叠矩阵
    pub fn smat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let smat = self.caler().smat(&form).clone();
        smat.into_pyarray(py).unbind()
    }

    /// 计算动能矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的动能矩阵
    pub fn tmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let tmat = self.caler().tmat(&form).clone();
        tmat.into_pyarray(py).unbind()
    }

    /// 计算势能（核吸引）矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个二维 NumPy 数组表示的势能矩阵
    pub fn vmat(&self, py: Python, form: String) -> Py<PyArray2<f64>> {
        let vmat = self.caler().vmat(&form).clone();
        vmat.into_pyarray(py).unbind()
    }

    /// 同时计算重叠、动能和势能三个单电子矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `cals`: 指定哪些矩阵需要计算，[calc_s, calc_t, calc_v]，1 表示计算，0 表示跳过
    /// - `form`: 积分格式，可选 "sph" 或 "car"
    ///
    /// ## 返回
    /// 返回一个包含三个二维 NumPy 数组的数组 [S, T, V]
    pub fn mat1e(&self, py: Python, cals: [usize; 3], form: String) -> [Py<PyArray2<f64>>; 3] {
        let [smat, tmat, vmat] = self.caler().mat1e(cals, &form);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }

    /// 把笛卡尔基下的单电子积分矩阵变换到球谐基
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `car_mat`: 笛卡尔基下的积分矩阵（二维 NumPy 数组）
    ///
    /// ## 返回
    /// 返回球谐基下的二维 NumPy 数组
    pub fn mat1e_sph(&self, py: Python, car_mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
        let car_mat = car_mat.readonly().as_array().to_owned();
        let mat = self.caler().mat1e_sph(&car_mat);
        mat.into_pyarray(py).unbind()
    }

    /// 把积分矩阵变换回未归一化的原始基
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `car_mat`: 笛卡尔基下的积分矩阵（二维 NumPy 数组）
    ///
    /// ## 返回
    /// 返回原始基下的二维 NumPy 数组
    pub fn mat1e_raw(&self, py: Python, car_mat: Bound<'_, PyArray2<f64>>) -> Py<PyArray2<f64>> {
        let car_mat = car_mat.readonly().as_array().to_owned();
        let mat = self.caler().mat1e_raw(&car_mat);
        mat.into_pyarray(py).unbind()
    }

    /// 在局域坐标系下同时计算三个单电子矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `stms`: 原子到壳层映射的字典（键为原子索引，值为 Stm 对象）
    /// - `cals`: 指定哪些矩阵需要计算，[calc_s, calc_t, calc_v]
    ///
    /// ## 返回
    /// 返回一个包含三个二维 NumPy 数组的数组 [S, T, V]
    pub fn mat1e_lcs(
        &self,
        py: Python,
        stms: HashMap<usize, Stm>,
        cals: [usize; 3],
    ) -> [Py<PyArray2<f64>>; 3] {
        let stms = stms.into_iter().map(|(idx, stm)| (idx, stm.core)).collect();
        let [smat, tmat, vmat] = self.caler().mat1e_lcs(&stms, cals);
        let smat = smat.into_pyarray(py).unbind();
        let tmat = tmat.into_pyarray(py).unbind();
        let vmat = vmat.into_pyarray(py).unbind();
        [smat, tmat, vmat]
    }
}

/// 双电子积分张量
///
/// 底层为 `rswfn::rsint::Mat2e`，按四个基函数下标 `(i, j, k, l)` 存取。
#[pyclass]
pub struct Mat2e {
    pub core: rswfn::rsint::Mat2e,
}

#[pymethods]
impl Mat2e {
    /// 创建一个双电子积分张量
    ///
    /// ## 参数
    /// - `nato`: 原子轨道（基函数）数量
    ///
    /// ## 返回
    /// 返回一个元素全为零的 Mat2e 对象
    #[new]
    pub fn new(nato: usize) -> Self {
        Self {
            core: rswfn::rsint::Mat2e::new(nato),
        }
    }

    /// 原子轨道（基函数）数量
    #[getter(nato)]
    pub fn get_nato(&self) -> usize {
        self.core.nato
    }

    /// 设置原子轨道（基函数）数量
    #[setter(nato)]
    pub fn set_nato(&mut self, nato: usize) {
        self.core.nato = nato;
    }

    /// 一维压缩存储的积分值
    #[getter(vals)]
    pub fn get_vals(&self) -> Vec<f64> {
        self.core.vals.clone()
    }

    /// 设置一维压缩存储的积分值
    #[setter(vals)]
    pub fn set_vals(&mut self, vals: Vec<f64>) {
        self.core.vals = vals;
    }

    /// 按四个基函数下标读取积分值
    ///
    /// ## 参数
    /// - `i`: 第一个基函数下标
    /// - `j`: 第二个基函数下标
    /// - `k`: 第三个基函数下标
    /// - `l`: 第四个基函数下标
    ///
    /// ## 返回
    /// 返回积分值 (ij|kl)
    pub fn get(&self, i: usize, j: usize, k: usize, l: usize) -> f64 {
        self.core.get(i, j, k, l)
    }

    /// 按四个基函数下标写入积分值
    ///
    /// ## 参数
    /// - `i`: 第一个基函数下标
    /// - `j`: 第二个基函数下标
    /// - `k`: 第三个基函数下标
    /// - `l`: 第四个基函数下标
    /// - `val`: 积分值 (ij|kl)
    pub fn set(&mut self, i: usize, j: usize, k: usize, l: usize, val: f64) {
        self.core.set(i, j, k, l, val);
    }

    /// 按壳层对下标读取积分值
    ///
    /// ## 参数
    /// - `p`: 第一个壳层对下标
    /// - `q`: 第二个壳层对下标
    ///
    /// ## 返回
    /// 返回该壳层对对应的积分值
    pub fn get_by_pair(&self, p: usize, q: usize) -> f64 {
        self.core.get_by_pair(p, q)
    }
}

/// 双电子积分引擎
///
/// 封装 `rswfn::rsint::Engine`，负责双电子积分的预筛与分块。
///
/// 由于上游 `Engine` 以借用方式持有 `Mole`，Python 对象中保存的是重建引擎
/// 所需的参数（`mole`、`form`、`screen`、`max_mem`），每次调用时按参数重建引擎，
/// 计算结果与一次构建完全一致。
#[pyclass]
pub struct Engine {
    pub mole: Mole,
    pub form: String,
    pub screen: f64,
    pub max_mem: usize,
    pub time_build: f64,
    pub nblock: usize,
    pub mem_bytes: usize,
}

impl Engine {
    /// 按当前参数重建底层积分引擎
    fn caler(&self) -> rswfn::rsint::Engine<'_> {
        rswfn::rsint::Engine::new(&self.mole.core, &self.form, self.screen, self.max_mem)
    }
}

#[pymethods]
impl Engine {
    /// 创建双电子积分引擎
    ///
    /// 创建时即构建一次引擎，并记录构建耗时、分块数、内存占用等信息。
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `form`: 基函数形式，可选 "sph" 或 "car"
    /// - `screen`: 积分预筛阈值
    /// - `max_mem`: 允许使用的最大内存（字节）
    ///
    /// ## 返回
    /// 返回一个 Engine 对象
    #[new]
    pub fn new(mole: Mole, form: String, screen: f64, max_mem: usize) -> Self {
        let (time_build, nblock, mem_bytes) = {
            let engine = rswfn::rsint::Engine::new(&mole.core, &form, screen, max_mem);
            (
                engine.time_build.as_secs_f64(),
                engine.nblock,
                engine.mem_bytes,
            )
        };
        Self {
            mole,
            form,
            screen,
            max_mem,
            time_build,
            nblock,
            mem_bytes,
        }
    }

    /// 引擎所属的分子对象
    #[getter(mole)]
    pub fn get_mole(&self) -> Mole {
        Mole {
            core: self.mole.core.clone(),
        }
    }

    /// 设置引擎所属的分子对象
    #[setter(mole)]
    pub fn set_mole(&mut self, mole: PyRef<'_, Mole>) {
        self.mole = Mole {
            core: mole.core.clone(),
        };
    }

    /// 基函数形式，可选 "sph" 或 "car"
    #[getter(form)]
    pub fn get_form(&self) -> String {
        self.form.clone()
    }

    /// 设置基函数形式
    #[setter(form)]
    pub fn set_form(&mut self, form: String) {
        self.form = form;
    }

    /// 积分预筛阈值
    #[getter(screen)]
    pub fn get_screen(&self) -> f64 {
        self.screen
    }

    /// 设置积分预筛阈值
    #[setter(screen)]
    pub fn set_screen(&mut self, screen: f64) {
        self.screen = screen;
    }

    /// 允许使用的最大内存（字节）
    #[getter(max_mem)]
    pub fn get_max_mem(&self) -> usize {
        self.max_mem
    }

    /// 设置允许使用的最大内存（字节）
    #[setter(max_mem)]
    pub fn set_max_mem(&mut self, max_mem: usize) {
        self.max_mem = max_mem;
    }

    /// 引擎构建耗时（秒）
    #[getter(time_build)]
    pub fn get_time_build(&self) -> f64 {
        self.time_build
    }

    /// 积分分块数量
    #[getter(nblock)]
    pub fn get_nblock(&self) -> usize {
        self.nblock
    }

    /// 引擎预估的内存占用（字节）
    #[getter(mem_bytes)]
    pub fn get_mem_bytes(&self) -> usize {
        self.mem_bytes
    }

    /// 检查积分分块信息
    ///
    /// ## 返回
    /// 返回 (分块数量, 原子轨道数量, 平均每个分块的积分数量) 三元组
    pub fn debug_check_blocks(&self) -> (usize, usize, f64) {
        self.caler().debug_check_blocks()
    }

    /// 原子轨道（基函数）数量
    ///
    /// ## 返回
    /// 返回基函数个数
    pub fn nato(&self) -> usize {
        self.caler().nato()
    }

    /// 是否为直接（direct）积分算法
    ///
    /// ## 返回
    /// 直接算法返回 True，否则返回 False
    pub fn is_direct(&self) -> bool {
        self.caler().is_direct()
    }

    /// 在局域坐标系下构建双电子积分张量
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `tmats`: 各原子的局域坐标变换矩阵列表
    /// - `pi_comp`: 每个原子的壳层对索引列表
    /// - `npi`: 壳层对数量
    ///
    /// ## 返回
    /// 返回构建好的 Mat2e 对象
    pub fn mat2e_local_pi(
        &self,
        py: Python,
        tmats: Vec<PyReadonlyArray2<f64>>,
        pi_comp: Vec<Vec<(usize, usize)>>,
        npi: usize,
    ) -> Mat2e {
        let tmats: Vec<Array2<f64>> = tmats.iter().map(|t| t.as_array().to_owned()).collect();
        let mole = self.mole.core.clone();
        let form = self.form.clone();
        let screen = self.screen;
        let max_mem = self.max_mem;
        let core = py.detach(move || {
            let engine = rswfn::rsint::Engine::new(&mole, &form, screen, max_mem);
            engine.mat2e_local_pi(&tmats, &pi_comp, npi)
        });
        Mat2e { core }
    }

    /// 计算库仑（J）与交换（K）矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `dmat`: 密度矩阵（二维 NumPy 数组）
    ///
    /// ## 返回
    /// 返回 (J 矩阵, K 矩阵) 二元组，均为二维 NumPy 数组
    pub fn jk(
        &self,
        py: Python,
        dmat: Bound<'_, PyArray2<f64>>,
    ) -> (Py<PyArray2<f64>>, Py<PyArray2<f64>>) {
        let dmat = dmat.readonly().as_array().to_owned();
        let mole = self.mole.core.clone();
        let form = self.form.clone();
        let screen = self.screen;
        let max_mem = self.max_mem;
        let (jmat, kmat) = py.detach(move || {
            let engine = rswfn::rsint::Engine::new(&mole, &form, screen, max_mem);
            engine.jk(&dmat)
        });
        (
            jmat.into_pyarray(py).unbind(),
            kmat.into_pyarray(py).unbind(),
        )
    }

    /// 计算开壳层的库仑（J）与自旋分辨交换（K）矩阵
    ///
    /// ## 参数
    /// - `py`: Python 解释器标记
    /// - `d_tot`: 总密度矩阵
    /// - `d_a`: alpha 自旋密度矩阵
    /// - `d_b`: beta 自旋密度矩阵
    ///
    /// ## 返回
    /// 返回 (J, K_alpha, K_beta) 三元组，均为二维 NumPy 数组
    pub fn jk_spin(
        &self,
        py: Python,
        d_tot: Bound<'_, PyArray2<f64>>,
        d_a: Bound<'_, PyArray2<f64>>,
        d_b: Bound<'_, PyArray2<f64>>,
    ) -> (
        Py<PyArray2<f64>>,
        Py<PyArray2<f64>>,
        Py<PyArray2<f64>>,
    ) {
        let d_tot = d_tot.readonly().as_array().to_owned();
        let d_a = d_a.readonly().as_array().to_owned();
        let d_b = d_b.readonly().as_array().to_owned();
        let mole = self.mole.core.clone();
        let form = self.form.clone();
        let screen = self.screen;
        let max_mem = self.max_mem;
        let (jmat, kamat, kbmat) = py.detach(move || {
            let engine = rswfn::rsint::Engine::new(&mole, &form, screen, max_mem);
            engine.jk_spin(&d_tot, &d_a, &d_b)
        });
        (
            jmat.into_pyarray(py).unbind(),
            kamat.into_pyarray(py).unbind(),
            kbmat.into_pyarray(py).unbind(),
        )
    }
}

/// 静电势（ESP）计算器
///
/// 封装 `rswfn::rsint::EspCalculator`。上游对象不持有 `Mole` 的所有权，
/// 因此 Python 对象中保存分子与密度矩阵，每次调用时按需重建计算器。
#[pyclass]
pub struct EspCalculator {
    pub mole: Mole,
    pub dmat: Array2<f64>,
}

impl EspCalculator {
    /// 按当前分子与密度矩阵重建底层计算器
    fn caler(&self) -> rswfn::rsint::EspCalculator {
        rswfn::rsint::EspCalculator::new(&self.mole.core, &self.dmat)
    }
}

#[pymethods]
impl EspCalculator {
    /// 创建静电势计算器
    ///
    /// ## 参数
    /// - `mole`: 分子对象
    /// - `dmat`: 密度矩阵（二维 NumPy 数组）
    ///
    /// ## 返回
    /// 返回一个 EspCalculator 对象
    #[new]
    pub fn new(mole: Mole, dmat: Bound<'_, PyArray2<f64>>) -> Self {
        let dmat = dmat.readonly().as_array().to_owned();
        Self { mole, dmat }
    }

    /// 静电势计算所用的分子对象
    #[getter(mole)]
    pub fn get_mole(&self) -> Mole {
        self.mole.clone()
    }

    /// 静电势计算所用的密度矩阵
    #[getter(dmat)]
    pub fn get_dmat(&self) -> Vec<Vec<f64>> {
        self.dmat.outer_iter().map(|row| row.to_vec()).collect()
    }

    /// 计算单个点的静电势
    ///
    /// ## 参数
    /// - `rc`: 点的笛卡尔坐标 [x, y, z]
    ///
    /// ## 返回
    /// 返回该点的静电势值
    pub fn val(&self, rc: [f64; 3]) -> f64 {
        self.caler().val(&rc)
    }

    /// 批量计算多个点的静电势
    ///
    /// ## 参数
    /// - `rcs`: 点坐标列表，每个元素为 [x, y, z]
    ///
    /// ## 返回
    /// 返回各点的静电势值列表
    pub fn vals(&self, rcs: Vec<[f64; 3]>) -> Vec<f64> {
        self.caler().vals(&rcs)
    }
}

/// 计算双电子积分在一维压缩存储中的下标
///
/// ## 参数
/// - `i`: 第一个基函数下标
/// - `j`: 第二个基函数下标
/// - `k`: 第三个基函数下标
/// - `l`: 第四个基函数下标
///
/// ## 返回
/// 返回对应的一维下标
#[pyfunction]
pub fn get_2eid(i: usize, j: usize, k: usize, l: usize) -> usize {
    rswfn::rsint::get_2eid(i, j, k, l)
}

/// 由积分张量与密度矩阵计算库仑（J）与交换（K）矩阵
///
/// ## 参数
/// - `py`: Python 解释器标记
/// - `dmat`: 密度矩阵（二维 NumPy 数组）
/// - `emat`: 双电子积分张量 Mat2e 对象
///
/// ## 返回
/// 返回 (J 矩阵, K 矩阵) 二元组，均为二维 NumPy 数组
#[pyfunction]
pub fn jk_mat(
    py: Python,
    dmat: Bound<'_, PyArray2<f64>>,
    emat: PyRef<'_, Mat2e>,
) -> (Py<PyArray2<f64>>, Py<PyArray2<f64>>) {
    let dmat = dmat.readonly().as_array().to_owned();
    let (jmat, kmat) = rswfn::rsint::jk_mat(&dmat, &emat.core);
    (
        jmat.into_pyarray(py).unbind(),
        kmat.into_pyarray(py).unbind(),
    )
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "rsint")?; // 为当前的rs文件创建一个子模块
    m.add_class::<RsInt>()?;
    m.add_class::<Mat2e>()?;
    m.add_class::<Engine>()?;
    m.add_class::<EspCalculator>()?;
    m.add_function(wrap_pyfunction!(get_2eid, &m)?)?;
    m.add_function(wrap_pyfunction!(jk_mat, &m)?)?;
    parent_module.add_submodule(&m) // 将子模块添加到父模块中
}
