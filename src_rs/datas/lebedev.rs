use pyo3::prelude::*;
use rswfn;

/// 根据 code 类型生成球面上的点、坐标和权重
#[pyfunction]
pub fn gen_oh(code: isize, a: f64, b: f64, v: f64) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    rswfn::datas::lebedev::gen_oh(code, a, b, v)
}

/// 将参数数组解析为多个 gen_oh 的输入，并合并结果为完整的 Lebedev 网格
#[pyfunction]
pub fn ld_fun(paras: Vec<[f64; 4]>) -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld_fun(&paras)
}

/// 生成6点 Lebedev 网格
#[pyfunction]
pub fn ld0006() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0006()
}

/// 生成14点 Lebedev 网格
#[pyfunction]
pub fn ld0014() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0014()
}

/// 生成26点 Lebedev 网格
#[pyfunction]
pub fn ld0026() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0026()
}

/// 生成38点 Lebedev 网格
#[pyfunction]
pub fn ld0038() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0038()
}

/// 生成50点 Lebedev 网格
#[pyfunction]
pub fn ld0050() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0050()
}

/// 生成74点 Lebedev 网格
#[pyfunction]
pub fn ld0074() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0074()
}

/// 生成86点 Lebedev 网格
#[pyfunction]
pub fn ld0086() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0086()
}

/// 生成110点 Lebedev 网格
#[pyfunction]
pub fn ld0110() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0110()
}

/// 生成146点 Lebedev 网格
#[pyfunction]
pub fn ld0146() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0146()
}

/// 生成170点 Lebedev 网格
#[pyfunction]
pub fn ld0170() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0170()
}

/// 生成194点 Lebedev 网格
#[pyfunction]
pub fn ld0194() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0194()
}

/// 生成230点 Lebedev 网格
#[pyfunction]
pub fn ld0230() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0230()
}

/// 生成266点 Lebedev 网格
#[pyfunction]
pub fn ld0266() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0266()
}

/// 生成302点 Lebedev 网格
#[pyfunction]
pub fn ld0302() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0302()
}

/// 生成350点 Lebedev 网格
#[pyfunction]
pub fn ld0350() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0350()
}

/// 生成434点 Lebedev 网格
#[pyfunction]
pub fn ld0434() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0434()
}

/// 生成590点 Lebedev 网格
#[pyfunction]
pub fn ld0590() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0590()
}

/// 生成770点 Lebedev 网格
#[pyfunction]
pub fn ld0770() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0770()
}

/// 生成974点 Lebedev 网格
#[pyfunction]
pub fn ld0974() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld0974()
}

/// 生成1202点 Lebedev 网格
#[pyfunction]
pub fn ld1202() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld1202()
}

/// 生成1454点 Lebedev 网格
#[pyfunction]
pub fn ld1454() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld1454()
}

/// 生成1730点 Lebedev 网格
#[pyfunction]
pub fn ld1730() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld1730()
}

/// 生成2030点 Lebedev 网格
#[pyfunction]
pub fn ld2030() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld2030()
}

/// 生成2354点 Lebedev 网格
#[pyfunction]
pub fn ld2354() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld2354()
}

/// 生成2702点 Lebedev 网格
#[pyfunction]
pub fn ld2702() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld2702()
}

/// 生成3074点 Lebedev 网格
#[pyfunction]
pub fn ld3074() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld3074()
}

/// 生成3470点 Lebedev 网格
#[pyfunction]
pub fn ld3470() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld3470()
}

/// 生成3890点 Lebedev 网格
#[pyfunction]
pub fn ld3890() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld3890()
}

/// 生成4334点 Lebedev 网格
#[pyfunction]
pub fn ld4334() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld4334()
}

/// 生成4802点 Lebedev 网格
#[pyfunction]
pub fn ld4802() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld4802()
}

/// 生成5294点 Lebedev 网格
#[pyfunction]
pub fn ld5294() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld5294()
}

/// 生成5810点 Lebedev 网格
#[pyfunction]
pub fn ld5810() -> (Vec<[f64; 3]>, Vec<f64>) {
    rswfn::datas::lebedev::ld5810()
}

pub fn register_module(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    let m = PyModule::new(parent_module.py(), "lebedev")?;
    m.add_function(wrap_pyfunction!(gen_oh, &m)?)?;
    m.add_function(wrap_pyfunction!(ld_fun, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0006, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0014, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0026, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0038, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0050, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0074, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0086, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0110, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0146, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0170, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0194, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0230, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0266, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0302, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0350, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0434, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0590, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0770, &m)?)?;
    m.add_function(wrap_pyfunction!(ld0974, &m)?)?;
    m.add_function(wrap_pyfunction!(ld1202, &m)?)?;
    m.add_function(wrap_pyfunction!(ld1454, &m)?)?;
    m.add_function(wrap_pyfunction!(ld1730, &m)?)?;
    m.add_function(wrap_pyfunction!(ld2030, &m)?)?;
    m.add_function(wrap_pyfunction!(ld2354, &m)?)?;
    m.add_function(wrap_pyfunction!(ld2702, &m)?)?;
    m.add_function(wrap_pyfunction!(ld3074, &m)?)?;
    m.add_function(wrap_pyfunction!(ld3470, &m)?)?;
    m.add_function(wrap_pyfunction!(ld3890, &m)?)?;
    m.add_function(wrap_pyfunction!(ld4334, &m)?)?;
    m.add_function(wrap_pyfunction!(ld4802, &m)?)?;
    m.add_function(wrap_pyfunction!(ld5294, &m)?)?;
    m.add_function(wrap_pyfunction!(ld5810, &m)?)?;
    parent_module.add_submodule(&m)
}
