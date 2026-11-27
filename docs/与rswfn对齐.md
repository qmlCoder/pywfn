# 与 rswfn 对齐说明

> 面向维护者的对齐手册。回答三个问题：**模块怎么对应**、**哪些差异是有意的**、**怎么验证还没对完**，最后附本轮新增接口清单。

| 项目 | 值 |
| --- | --- |
| 上游依赖 | `rswfn = "=0.2.0"`（`Cargo.toml`） |
| pywfn 版本 | `1.0.22`（`pyproject.toml`） |
| 扩展模块 | `pywfn._core`（`module-name = "pywfn._core"`，`maturin` + `pyo3/abi3-py310`） |
| 文档基准 | rswfn rustdoc 快照（`target/apialign/rswfn_api.json`） |
| 本文快照时间 | 2026-11-27 |

## 0. 三层结构与两条约定

```
rswfn (Rust crate, 0.2.0)
  └── src_rs/*.rs            pyo3 绑定层：只做类型转换与转发
        └── _core.pyd        pywfn._core：maturin 构建出的扩展模块
              └── src_py/pywfn/*.py   最终使用的 Python 包：类型标注 + 中文文档字符串
```

两条必须遵守的约定（来自 `AGENTS.md`）：

1. **`src_rs/` 只翻译**：不含算法、不写注释，所有计算都回到 `rswfn`。
2. **`src_py/pywfn/` 只封装**：每个封装类持有底层对象 `self.core`，每个方法转调它；参数与返回值都做 NumPy / Python 类型标注。

一个容易踩坑的细节：`src_py/pywfn/_core.pyi` 是 **空文件**（0 字节），所以 `_core` 对类型检查器完全不可见，`src_py` 里大量 `# type: ignore` 是在补这个洞，不是历史遗留垃圾。

## 1. 总体结构

### 1.1 模块对应关系

| rswfn 模块 | pywfn Python 模块 | `_core` 子模块 | 说明 |
| --- | --- | --- | --- |
| `base` | `pywfn.base`（`atoms` / `basis` / `coefs` / `mole` / `vect`） | `_core.base` | 12 个结构体：`Ato` `Atom` `Atoms` `Atos` `Bas` `Basis` `Coefs` `Mole` `Shl` `Shls` `Stm` `Vec3` |
| `datas` | `pywfn.datas` | `_core.datas` → `btrans` / `consts` / `dhmo2` / `elems` / `integ` / `lebedev` / `rdens` / `roots` / `rys_root` / `temps` | 包级再导出 `FACS` 与 `Rhm` / `Whm`；`datas/march.py` 是历史遗留的 `marchData` 查表数据，与 `maths.march` 的算法无关 |
| `editor` | `pywfn.editor` | `_core.editor` | `Editor` 类 + `rotate_points` |
| `loger` | `pywfn.loger` | `_core.loger` | `init_log` / `init_stderr_log` |
| `maths` | `pywfn.maths` | `_core.maths` → `march`（`mc` / `mt` / `sn`）、`matrix`、`optim`（`cg` / `dn` / `nt` / `qn` / `sd`） | 顶层 21 个函数 |
| `props.atomprop` | `pywfn.props.atomprop` | `_core.props.atomprop` → `activity` / `charge` / `direction` / `resp` | |
| `props.bondprop` | `pywfn.props.bondprop` | `_core.props.bondprop` → `direction` / `energy` / `order` | Python 侧另有 `selectPI.py` |
| `props.fragprop` | `pywfn.props.fragprop` | `_core.props.fragprop` → `energy` | 目前只有 `Calculator` |
| `props.gridprop` | `pywfn.props.gridprop` | `_core.props.gridprop` → `density` / `dftgrid` / `potential` / `topolog` / `wfnfunc` | 顶层四个网格类；Python 侧另有 `hmowfn.py` |
| `props.moleprop` | `pywfn.props.moleprop` | `_core.props.moleprop` → `aromacity` / `energy` / `newmole` / `population` / `variants` | `newmole` 是旧名，`variants` 与 rswfn 对齐，两者目前都保留 |
| `props.orbtprop` | `pywfn.props.orbtprop` | `_core.props.orbtprop` → `mocv` / `obtmat` / `obtval` | 原来的 `datas/mocv` 已迁到这里（见 4.5） |
| `qmcal` | `pywfn.qmcal` | `_core.qmcal` → `dhmo` / `dhmo2` / `diis` / `guess` / `hf` / `ks` / `lks` / `se` | 模块级另有 `SlaterX`；`ks.B3LYP_EXACT_X` 为 float 常量 |
| `reader` | `pywfn.reader` | `_core.reader` | 13 个具体 Reader + `FileType` 等；`trait Reader` 不导出（2.1） |
| `rsbse` | `pywfn.rsbse` | `_core.rsbse` | 基组下载/解析；`BASIS_SETS` 以 list 形式暴露 |
| `rsint` | `pywfn.rsint` | `_core.rsint` | `RsInt` / `Engine` / `EspCalculator` / `Mat2e` + `get_2eid` / `jk_mat` |
| `rsxcs` | `pywfn.rsxcs` | `_core.rsxcs` | `XcComponent` / `XcFunc` / `RsXcs` / `RksEval` / `UksEval` + 4 个函数 |
| `space` | `pywfn.space` | `_core.space` | `a2m_weits` / `calc_EIEBA` / `get_grids` |
| `tools` | `pywfn.tools` | `_core.tools` → `log` | `_core` 只有 `log.Tool`；Python 侧的 `gjf.py` / `xml.py` / `xmlPES.py` 是 pywfn 自有扩展，rswfn 文档中没有同名模块 |
| `utils` | `pywfn.utils` | `_core.utils` | 11 个函数，其中 `rswfn_version` 是 pywfn 额外加的 |
| `writer` | `pywfn.writer` | `_core.writer` | 7 个具体 Writer；`trait Writer` 不导出（2.1） |

### 1.2 注册与导出入口

| 层 | 位置 | 行为 |
| --- | --- | --- |
| Rust | `src_rs/lib.rs` | `#[pymodule] fn _core(...)` 依次注册 15 个子模块：`base` `datas` `editor` `loger` `maths` `props` `qmcal` `reader` `rsbse` `rsint` `rsxcs` `space` `tools` `utils` `writer` |
| Python | `src_py/pywfn/__init__.py` | 再导出其中 14 个（不含 `tools`），`__all__` 同名列；另有 `rswfn_version()` 与 `__rswfn_version__` |

对应关系：`pywfn.<一级>.<二级>.<名字>` ↔ rswfn 的 `<一级>::<二级>::<名字>`。`_core` 是这条路径的"平铺版本"，差异见 2.4。

## 2. 有意保留的差异

下面五处差异**不是遗漏**，改动前请先确认真的需要改。

### 2.1 `trait Reader` / `trait Writer` 只导出具体类型

rswfn 的 `reader::Reader`、`writer::Writer` 是 trait（`target/apialign/rswfn_api.json` 中 `kind = "trait"`）。pyo3 只能导出具体类型：trait 不能写成 `#[pyclass]`，也无法在 Python 侧充当抽象基类。

因此：

- `_core.reader` 暴露 13 个具体 Reader（`CubReader` `DmoReader` `FchReader` `GjfReader` `InpReader` `LogReader` `MdeReader` `MfnReader` `MolReader` `NonReader` `OutReader` `SdfReader` `XyzReader`）以及 `get_fileType` / `charge_number` / `read_lines` / `shl_basis` / `shl_forms` / `parse_job`；
- `_core.writer` 暴露 7 个具体 Writer（`CubWriter` `DmoWriter` `FchWriter` `GjfWriter` `MfnWriter` `ObjWriter` `XyzWriter`）；
- Python 侧靠**鸭子类型**统一：每个 Reader 都实现 `get_atoms()` / `get_basis()` / `get_coefs()`，`pywfn.reader.get_reader(path)` 按后缀分派（映射表 `_READERS`，见 `src_py/pywfn/reader/__init__.py`）。

同理的还有 `reader::BasData`：rswfn 里是类型别名 `(usize, usize, usize, f64, f64)`，`_core` 里变成带 `as_tuple()` 的类。

### 2.2 缓存句柄静态量不导出，改为"三件套"

rswfn 的 `datas` 中用 `static`（锁保护的映射）缓存高频查询的表。这类句柄跨 FFI 暴露会牵扯锁与 GIL 语义，且属于 Rust 侧实现细节，因此**没有**出现在 `_core` 里。取而代之的是同名函数的三个版本：

| rswfn static | 缓存版 | 不走缓存 | 预热缓存 |
| --- | --- | --- | --- |
| `btrans::GET_CAR_SYMS` | `get_car_syms` | `get_car_syms_no_cache` | `get_car_syms_prime_cache` |
| `btrans::GET_MIX_SYMS` | `get_mix_syms` | `get_mix_syms_no_cache` | `get_mix_syms_prime_cache` |
| `btrans::GET_SPH_SYMS` | `get_sph_syms` | `get_sph_syms_no_cache` | `get_sph_syms_prime_cache` |
| `btrans::LEN2SYMS` | `len2syms` | `len2syms_no_cache` | `len2syms_prime_cache` |
| `btrans::SYM2ANG` | `sym2ang` | `sym2ang_no_cache` | `sym2ang_prime_cache` |
| `btrans::SYM2LMN` | `sym2lmn` | `sym2lmn_no_cache` | `sym2lmn_prime_cache` |
| `elems::ATM2RAD` | `atm2rad` | `atm2rad_no_cache` | `atm2rad_prime_cache` |
| `elems::ATM2ROW` | `atm2row` | `atm2row_no_cache` | `atm2row_prime_cache` |
| `elems::ATM2SYM` | `atm2sym` | `atm2sym_no_cache` | `atm2sym_prime_cache` |
| `elems::GET_BY_IDX` | `get_by_idx` | `get_by_idx_no_cache` | `get_by_idx_prime_cache` |
| `elems::GET_BY_SYM` | `get_by_sym` | `get_by_sym_no_cache` | `get_by_sym_prime_cache` |
| `elems::ZEFF` | `zeff` | `zeff_no_cache` | `zeff_prime_cache` |

共 12 组、36 个函数，全部已暴露。核查方式：

```python
from pywfn import _core
assert not hasattr(_core.datas.btrans, "GET_CAR_SYMS")   # 缓存句柄：不导出
assert hasattr(_core.datas.btrans, "get_car_syms_no_cache")  # 三个函数：导出
```

注意区分 **缓存句柄** 与 **只读数据表**：后者照常暴露。`DMAT` / `FMAT` / `GMAT` / `HMAT`、`RHM` / `WHM`、`GJF` / `SI` / `XYZ`、`ADP` / `RAD_DENS` / `RAD_DIST`、`MARCH`、`ELEMS` / `ELEM_DATA` / `ELES` / `NLMS` / `OXIS` 在 `_core` 中都是**函数**（调用一次返回一份数据，避免外部改到内部表）；个别量直接给值：`_core.qmcal.ks.B3LYP_EXACT_X`（`float`）、`_core.rsbse.BASIS_SETS`（`list`）。

### 2.3 需要闭包参数的 API：改收 Python 可调用对象

| rswfn | pywfn 签名 |
| --- | --- |
| `maths::march::extract_region` | `extract_region(func, cube, isov, step, xyzs=None, limit=None, extract=None)`；`cube` 收 Python 的 `Cube` 对象 |
| `maths::march::mc::march_cube_v2` | `march_cube_v2(func, cube, isov, step, xyzs=None, limit=None)` |
| `maths::march::mt::march_tetra_v2` | `march_tetra_v2(func, cube, isov, step, xyzs=None, limit=None)` |
| `maths::march::smooth_isosurface` | `smooth_isosurface(verts, types, isov, func, iterations)`（回调在最后几参） |
| `maths::march::snap_to_isosurface` | `snap_to_isosurface(verts, types, isov, func, h)` |
| `maths::optim::{cg,dn,nt,qn,sd}::minimize` | `minimize(func, paras: list[float], conv_thresh: float)` |

实现方式：`Py<PyAny>` 在 Rust 里包成闭包（`src_rs/maths/march.rs` 的 `wrap_func`、`src_rs/maths/optim/mod.rs` 的 `parse_funceval`），求值时用 `Python::attach` 重新拿 GIL 调回 Python。

两点注意事项：

- 回调**每次求值都要重新取 GIL**，`extract_region` / `minimize` 内部的并行加速会被回调本身吃掉；能用 NumPy 向量化预计算的，不要放在回调里逐点算。
- `minimize` 的回调返回值支持三种写法：直接返回 `float`、返回 `FuncEval`、返回 `(value, grad, hessian)` 元组。这样只给函数值时也能跑（走数值差分），给了梯度/海森时精度更好。

Python 侧在 `pywfn.maths.optim` 上再包了一层：`minimize(func, paras, conv_thresh=1e-6, method="cg")` 用 `method` 字符串在 `cg` / `dn` / `nt` / `qn` / `sd` 之间选。

### 2.4 `_core` 内部平铺，Python 包层再分层

`_core` 保留 rswfn 的**模块 → 子模块**两级（`_core.datas.btrans`、`_core.maths.march.mc`、`_core.props.atomprop.charge`），但凡是 rswfn 里"一个文件一个类型"的 `reader` / `writer`，在 `_core` 中被拍平成一层：

| rswfn 路径 | `_core` 路径 | pywfn Python 路径 |
| --- | --- | --- |
| `reader::fch::FchReader` | `_core.reader.FchReader` | `pywfn.reader.fch.FchReader` |
| `writer::fch::FchWriter` | `_core.writer.FchWriter` | `pywfn.writer.fch.FchWriter` |
| `maths::march::mc::march_cube_v1` | `_core.maths.march.mc.march_cube_v1` | `pywfn.maths.march.march_cube_v1` |
| `props::orbtprop::mocv::MocvDest` | `_core.props.orbtprop.mocv.MocvDest` | `pywfn.props.orbtprop.mocv.MocvDest` |

Python 包层负责把文件层级补回来：一个类一个同名 `.py`，`__init__.py` 里做汇总导出。所以**路径不要按字符串直接比对**，要按 1.1 的映射表换算（`diff_api.py` 的 `RENAMES` 就是干这个的）。

### 2.5 非 `Clone` / 带生命周期的类型：存字段、调用时重建

`rswfn` 有些类型不是 `Clone`，或带 `'_` 生命周期，pyo3 的 `#[pyclass]` 不能直接持有它们。做法是：**只保存重建所需的字段**，每次调用时临时构造底层对象。

| pywfn 类 | 保存的字段 | 重建方式 | 代价 |
| --- | --- | --- | --- |
| `maths.march.Cube` | `p0`、`p1` | `rswfn::maths::march::Cube::new(p0, p1)` | 可忽略 |
| `rsint.RsInt` | `mole` | `rswfn::rsint::RsInt<'_>::new(&mole.core)` | 中：`smat` / `tmat` / `vmat` / `mat1e*` 每次都重建 |
| `rsint.Engine` | `mole`、`form`、`screen`、`max_mem`（另存构造时的 `time_build` / `nblock` / `mem_bytes` 供属性读取） | `rswfn::rsint::Engine<'_>::new(&mole.core, &form, screen, max_mem)` | **高**：除 `time_build` / `nblock` / `mem_bytes` 三个属性外，每个方法（`debug_check_blocks` / `nato` / `is_direct` / `mat2e_local_pi` / `jk` / `jk_spin`）都会重建积分引擎；不要逐点调用，攒成一批 |
| `rsint.EspCalculator` | `mole`、`dmat` | `rswfn::rsint::EspCalculator::new(&mole.core, &dmat)` | 中：每次 `val` / `vals` 重建 |
| `props.orbtprop.mocv.MocvDest` | `kind`（Python 侧的镜像枚举 `MocvKind`） | `to_core()` 转成 rswfn 的 `MocvDest` 枚举 | 可忽略 |

这类类的 Python 封装遵循同一模式：`__init__` 里构造一次（或调用时重建），属性读取从保存的字段走。新增同类包装时请照此办理，并在类的 docstring 里标明"每次调用会重建底层对象"。

## 3. 如何校验对齐

思路：**rswfn 的 rustdoc 是唯一基准**（按项目约定，对齐时不读 rswfn 源码，只看文档接口），把文档抽成 JSON，再和 `_core` 的运行时内省结果做集合差。

### 3.1 前置：先有 rustdoc

在 rswfn 仓库执行 `cargo doc --no-deps`，产物默认在 `E:\code\rswfn\target\doc\rswfn`。路径不同时用环境变量覆盖：

```powershell
$env:RSWFN_DOC = "D:\code\rswfn\target\doc\rswfn"
```

### 3.2 `target/apialign/` 下的脚本

这些脚本**存在于本机**（`target/` 在 `.gitignore` 中，不进版本库，`cargo clean` 或换机器会丢）。其中 `diff_api.py` / `diff_methods.py` 要内省当前已构建的 `_core`，所以先 `maturin develop`；`extract_api.py` 只读 rustdoc，不依赖构建。

| 脚本 | 作用 | 用法 |
| --- | --- | --- |
| `extract_api.py` | 扫 rustdoc HTML，抽出模块 / 结构体 / 枚举 / 函数 / 方法签名 → `rswfn_api.json` | `python target/apialign/extract_api.py`（可用 `RSWFN_DOC` / `RSWFN_OUT` 覆盖路径） |
| `show_api.py` | 按前缀打印接口规格（声明 + 方法 + 变体），写绑定时照抄 | `python target/apialign/show_api.py props.orbtprop` |
| `diff_api.py` | `_core` 内省 vs `rswfn_api.json`，输出「缺失模块 / 缺失条目 / pywfn 多出条目」 | `python target/apialign/diff_api.py` |
| `diff_methods.py` | 逐结构体比对**方法名**，只报 pywfn 还缺哪些 | `python target/apialign/diff_methods.py` |
| `diff_docs.py` | 文档版 rswfn 与**实际依赖版本**（`dep_api.json`）的声明/方法差异，用来确认"文档里有、能绑定的也有" | `python target/apialign/diff_docs.py` |
| `fix_insert.py` | 一次性修补脚本（把误插进函数体的方法块搬回 `impl`），正常流程用不到 | 不建议再跑 |
| `spec.md` | `show_api.py` 的完整输出快照（约 69 KB），可当离线规格书读 | 直接打开 |

### 3.3 推荐流程

1. `maturin develop` 构建当前 `_core`；
2. `python target/apialign/extract_api.py` 刷新 `rswfn_api.json`；
3. `python target/apialign/diff_api.py` → 先消灭「缺失模块」，再逐条看「缺失条目」；
4. `python target/apialign/diff_methods.py` → 补结构体方法；
5. 抽查几条 `show_api.py <前缀>` 对签名（参数名、可选参数、返回类型）；
6. `pytest tests/test_align.py -q` 跑冒烟测试。

### 3.4 用 Python 内省做定点核查

`diff_api.py` 的核心就是下面这段，临时查单个模块时可以直接抄：

```python
import types
from pywfn import _core

def walk(mod, prefix=""):
    out = {}
    for name in sorted(dir(mod)):
        if name.startswith("__"):
            continue
        obj = getattr(mod, name)
        full = f"{prefix}{name}"
        if isinstance(obj, types.ModuleType):
            out.update(walk(obj, full + "."))
        else:
            out[full] = obj
    return out

names = walk(_core)
print([n for n in names if n.startswith("props.orbtprop.mocv")])
```

查方法与签名：

```python
import inspect
from pywfn import _core
print([m for m in dir(_core.rsint.Engine) if not m.startswith("_")])
print(inspect.signature(_core.maths.march.extract_region))
```

### 3.5 判读规则：哪些差异是预期的

看到 `diff_api.py` 报缺失时，先按这张表过滤，避免白干：

| 差异 | 是否预期 | 处理 |
| --- | --- | --- |
| `reader::Reader`、`writer::Writer` | 预期 | 不导出（2.1） |
| `GET_CAR_SYMS` / `ATM2RAD` / `ZEFF` / `LEN2SYMS` 等 static | 预期 | 已换成三件套（2.2），用 `*_no_cache` / `*_prime_cache` 对账 |
| `reader::FchReader` 与 `_core.reader.FchReader` 路径不等 | 预期 | 平铺（2.4），按映射表比对 |
| 用闭包参数的函数（`extract_region`、`minimize` …）签名不同 | 预期 | Python 侧收可调用对象（2.3） |
| 非 `Clone` 类型的字段/方法看起来"少了" | 先查 2.5 | 可能只保留了重建字段 |
| 其他 | **真缺口** | 在 `src_rs` 补绑定 → `src_py` 补封装 → 补测试 |

### 3.6 注意事项

- `diff_api.py` 里的 `RENAMES` 记录"同名不同路径"的重命名，改路径时要同步维护，目前是：`datas.rsroot → datas.rys_root`、`datas.mocv → props.orbtprop.mocv`、`props.moleprop.newmole → props.moleprop.variants`、`props.gridprop.potential → props.gridprop.dftgrid`。
- 脚本里写的是绝对路径 `E:\code\pywfn\target\apialign\rswfn_api.json`，换工作目录要改。
- `rswfn_api.json` 与 `dep_api.json` 是**快照**：升级 rswfn 后必须重新生成，否则 `diff_docs.py` 比的是旧文档。
- `diff_api.py` 只看名字层面，不校验参数类型；签名以 `show_api.py` 和 rustdoc 为准。

## 4. 新增/变更的接口清单

统计口径：相对上次提交 `HEAD`。本轮 `src_py/` 新增 22 个文件、修改 35 个；`src_rs/` 新增 26 项（含 `rsbse/` `rsxcs/` `maths/optim/` 三个目录）、修改 49 个、删除 2 个（`datas/mocv.rs`、`datas/rsroot.rs`）。

### 4.1 新增的 Python 模块文件（22 个）

| 文件 | 主要内容 |
| --- | --- |
| `src_py/pywfn/loger.py` | `init_log` / `init_stderr_log` |
| `src_py/pywfn/space.py` | `a2m_weits` / `calc_EIEBA` / `get_grids` |
| `src_py/pywfn/utils.py` | 11 个工具函数（`file_name` `file_type` `fort_form` `all_same` `symm_matrix` `posi_matrix` `has_nan` `same_matrix` `show_matrix` `comp_matrix` `rswfn_version`） |
| `src_py/pywfn/rsbse.py` | `BASIS_SETS`、`from_file` / `from_name` / `from_text` / `read_basis` / `supported_names` |
| `src_py/pywfn/rsxcs.py` | `XcComponent` / `XcFunc` / `RsXcs` / `RksEval` / `UksEval`、`eval_rks` / `eval_uks` / `supported` / `supported_names` |
| `datas/dhmo2.py` | `get_paras` |
| `datas/integ.py` | `RHM` / `WHM`（模块级 `ndarray`） |
| `datas/rdens.py` | `get_rad_dens_v1` / `get_rad_dens_v2`（`ADP` / `RAD_DENS` / `RAD_DIST` 以函数暴露于 `_core`） |
| `datas/roots.py` | `hermite` / `legendre` / `dub_roots` / `rys_roots` + `RYS_ALPS`、`RYS_ROOTS_1..9`、`RYS_WEITS_1..9` |
| `datas/rys_root.py` | `rys_root1` ~ `rys_root5`、`rys_roots` + `POLY_*` 六张表 |
| `maths/optim.py` | `FuncEval`、`minimize`、`cg` / `dn` / `nt` / `qn` / `sd`、`OptCallable` |
| `props/atomprop/resp.py` | `Calculator`（`esp_values` `mk_points` `resp` `resp_with`）、`EspType`、`RespOptions`、`RespResult` |
| `props/bondprop/energy.py` | `BondEnergy`、`Calculator`（`bond_energy` `bond_parts`） |
| `props/gridprop/topolog.py` | `Calculator`（`analysis` `get_func` `get_grids` `get_point`） |
| `props/moleprop/variants.py` | `Calculator`（`mole_spin` `mole_pi_mocv` `mole_pi_pocv`） |
| `props/orbtprop/mocv.py` | `OrderDType`、`MocvDest`、`mocv` / `mocv_atos` / `mocv_syms`、`get_nmat` / `get_tmat` / `get_tmat_sph` / `get_tmat_car`、`tmat_s` / `tmat_p` / `tmat_d` / `tmat_f` / `tmat_g`、`nmat_d` / `nmat_f` / `nmat_g` |
| `qmcal/dhmo.py` | `DmoRes`、`dhmo_calc` / `effective_charge` / `free_valence` / `output` / `overlap_integral` / `HMO_Analysis` |
| `qmcal/dhmo2.py` | `dhmo_calc_v2` / `get_integ_d` / `get_integ_x` |
| `qmcal/diis.py` | `Diis`（`commutator` `extrapolate` `reset`） |
| `qmcal/guess.py` | `ao_owners` / `sad_dmat` / `use_sad_enabled` |
| `qmcal/ks.py` | `RKS` / `UKS` / `RKSResult` / `UKSResult` / `GridBlock` / `SpinRho`、`rks_ele_eng` / `uks_ele_eng` |
| `qmcal/lks.py` | `LKS` / `LKSResult` |

### 4.2 包级导出变化

| 位置 | 变化 |
| --- | --- |
| `pywfn/__init__.py` | 新增 14 个子模块的导入与 `__all__`（`base` `datas` `editor` `loger` `maths` `props` `qmcal` `reader` `rsbse` `rsint` `rsxcs` `space` `utils` `writer`，不含 `tools`） |
| `src_rs/lib.rs` | 新增 `rsbse`、`rsxcs` 两个子模块的 `register_module` 调用 |
| `datas/__init__.py` | 新增 `FACS`（阶乘表，来自 `_core.datas.FACS()`）；`Rhm` / `Whm` 改成从 `datas/integ.py` 取值 —— 原来这两个 10×10 表是手写在 Python 里的，现在统一取自 `_core.datas.integ` |
| `props/gridprop/__init__.py` | 新增 `GridVals` 导出 |
| `qmcal/__init__.py` | 新增 `SlaterX` |
| `props/moleprop/energy.py` | 新增 `PiEnergy` |
| `base/__init__.py` | 导出新增 `Ato`、`Bas` |

### 4.3 已有模块新增的类 / 函数

| 模块 | 新增 |
| --- | --- |
| `pywfn.base` | `Ato`、`Bas` 两个封装类 |
| `pywfn.datas.btrans` | `lmn2sym` / `sym2lmn` / `sym2ang` / `sym2stp` / `ang2lmns` / `ang2syms` / `get_fmat`，加上 6 组 `*_no_cache` / `*_prime_cache` |
| `pywfn.datas.elems` | 6 组 `*_no_cache` / `*_prime_cache`（`atm2rad` `atm2row` `atm2sym` `get_by_idx` `get_by_sym` `zeff`） |
| `pywfn.datas.temps` | `GJF` / `SI` / `XYZ` |
| `pywfn.datas.consts` | `BOHR` 改为取自 `_core.datas.consts.BOHR()`（原来在 Python 里手写 `1.0/0.52917721092`） |
| `pywfn.maths` | 15 个顶层函数：`arg_max` `arg_min` `bolzman_dist` `calc_dist` `det` `eigh` `get_plan_by_3points` `inv` `sign` `vec_cross` `vec_dist` `vec_len` `vec_norm` `vec_project` `vec_psd` |
| `pywfn.maths.march` | `Cube`；`extract_region` `smooth_isosurface` `snap_to_isosurface` `march_cube_v2` `march_tetra_v1` `march_tetra_v2` `surface_nets_v1` `flatten` `grids2voxel` `voxel2verts` `merge_verts` `value_merge` `voxel2cubes` |
| `pywfn.props.atomprop` | `activity.element_valence`；`charge.hirshfeld` / `charge.qmat2eles`；`direction` 的 6 个方法；`resp` 整模块 |
| `pywfn.props.bondprop` | `order.order()`；`energy` 整模块 |
| `pywfn.props.gridprop` | `GridVals`；`dftgrid` 的 `atm_grids_meta` / `atm_grids_p` / `elm_grids` / `mol_grids_meta` / `mol_grids_p`；`density` 的 9 个方法（`activate` `curvature` `ele_esp` `tot_esp` `nuc_esp` `mol_rho_spin` `mole_spin` `mole_pi_mocv` `mole_pi_pocv`）；`wfnfunc.ato_wfn` / `wfnfunc.stf`；`topolog` 整模块 |
| `pywfn.props.moleprop` | `energy` 的 `PiEnergy` 与 4 个 `pi_*` 方法；`variants` 整模块 |
| `pywfn.props.orbtprop` | `mocv` 整模块；`obtmat.Calculator.LMO`；`obtval.Calculator.eles` |
| `pywfn.qmcal` | `hf.UHF.calc_udmat` / `calc_uhf`；`se.RSE.calc_rse`；6 个新模块（4.1） |
| `pywfn.reader` / `pywfn.writer` | Python 层本轮没有新增模块文件，增量主要在 `_core`（见 4.4） |

### 4.4 已有类补齐的方法（精选）

| 类 | 新增方法 |
| --- | --- |
| `base.Stm` | `blank` `check` `copy` `dirs` `from_mat` `rand_ex` `rand_ez` `relative` `vals` |
| `base.Atom` | `irow` `vec` |
| `base.Atoms` | `ang` `bonds` `cent` `dmat` `equal` `idxs` `nuc_eng` `rads` `set_bonds` `stms` `xyzs_ang` |
| `base.Basis` | `ato_keys` `build` `equal` `get_by_atm` `int_form` `is_none` `nbas` `nshl` `unify_form` |
| `base.Coefs` | `equal` `is_none` `obt_eles` |
| `base.Mole` | `build` `dist` `equal` `formula` `is_line` `ncrg` `rings` |
| `_core.reader.*` | 各 Reader 的 `read_*`（例如 `LogReader.read_geomes` / `read_freqs` / `read_nics` / `read_energys` / `read_smat` / `read_tmat` / `read_vmat`、`FchReader.read_crgs` / `read_stms`、`GjfReader.read_job` / `read_gau_job`）、`MfnReader`、`FileType`、`GauJob`、`Marks`、`Vibe`、`BasData`、`parse_job`、`charge_number`、`read_lines`、`shl_basis`、`shl_forms` |
| `_core.writer.*` | `MfnWriter`；各 Writer 的 `build` / `from_mole` / `read_mole` / `set` / `save`（具体有哪几个按类型而异，`ObjWriter` 只有 `save`） |
| `_core.rsint` | `Engine`、`EspCalculator`、`Mat2e`、`get_2eid`、`jk_mat`；`RsInt.mat1e_lcs` / `mat1e_raw` / `mat1e_sph` |
| `_core.maths.optim` | `FuncEval` 与 5 个 `minimize`（按算法分子模块） |

### 4.5 迁移与删除

| 变化 | 说明 |
| --- | --- |
| `datas/mocv`（Python + Rust）删除 | 轨道系数变换并入 `props.orbtprop.mocv`，与 rswfn 的 `props::orbtprop::mocv` 对齐 |
| `datas/rstoot.py` / `datas/rsroot.rs` 删除 | 重命名为 `datas/rys_root`，与 rswfn 的 `datas::rys_root` 对齐 |
| `props/moleprop/newmole` 保留 | 与 `variants` 内容重复但 `_core` 两套都在；新增代码请用 `variants`（对齐名） |

### 4.6 当前落地状态

| 项 | 状态 |
| --- | --- |
| pywfn.reader 的重导出 | 13 个 Reader 与 get_reader/supports，以及 _core.reader 的 FileType / BasData / Marks / Vibe / GauJob / charge_number / get_fileType / 
ead_lines / shl_basis / shl_forms / parse_job 全部转出 |
| pywfn.writer 的重导出 | 7 个 Writer（含 MfnWriter） |
| pywfn.rsint | RsInt / Mat2e / Engine / EspCalculator + get_2eid / jk_mat |
| _core.pyi | 仍为空文件，_core 无类型桩，故 src_py 中保留大量 # type: ignore |
| 测试 | pytest tests -q → 22 passed（含 tests/test_align.py 共 11 项） |

## 5. 维护动作清单

新增/对齐一个接口时的固定顺序：

1. 读 rustdoc（`show_api.py <前缀>`），确认签名与返回类型；
2. `src_rs/<模块>.rs` 加 `#[pyclass]` / `#[pyfunction]` / `#[pymethods]`，并在对应 `register_module` 里注册；
3. 若是新的一级模块，同时改 `src_rs/lib.rs` 的 `pub mod` 与注册调用；
4. `src_py/pywfn/<模块>/<文件>.py` 加封装（`self.core` 转发 + 类型标注 + 中文 docstring，注释从 rswfn 文档翻译）；
5. 新模块要加进 `pywfn/__init__.py` 的导入与 `__all__`；
6. `tests/` 补冒烟测试（不跑 SCF）；
7. `maturin develop` → `pytest -q` → `diff_api.py` / `diff_methods.py` 复查；
8. 发版流程按 `AGENTS.md`：改版本号 → `maturin build --release` → PyPI + GitHub release。
