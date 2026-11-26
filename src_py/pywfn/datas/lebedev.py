"""
记录 Lebedev 网格数据
"""

from pywfn import _core

core = _core.datas.lebedev  # type: ignore


def gen_oh(
    code: int, a: float, b: float, v: float
) -> tuple[list[float], list[float], list[float], list[float]]:
    """根据 code 类型生成球面上的点、坐标和权重"""
    return core.gen_oh(code, a, b, v)  # type: ignore


def ld_fun(paras: list[list[float]]) -> tuple[list[list[float]], list[float]]:
    """将参数数组解析为多个 gen_oh 的输入，并合并结果为完整的 Lebedev 网格"""
    return core.ld_fun(paras)  # type: ignore


def ld0006() -> tuple[list[list[float]], list[float]]:
    """生成6点 Lebedev 网格"""
    return core.ld0006()  # type: ignore


def ld0014() -> tuple[list[list[float]], list[float]]:
    """生成14点 Lebedev 网格"""
    return core.ld0014()  # type: ignore


def ld0026() -> tuple[list[list[float]], list[float]]:
    """生成26点 Lebedev 网格"""
    return core.ld0026()  # type: ignore


def ld0038() -> tuple[list[list[float]], list[float]]:
    """生成38点 Lebedev 网格"""
    return core.ld0038()  # type: ignore


def ld0050() -> tuple[list[list[float]], list[float]]:
    """生成50点 Lebedev 网格"""
    return core.ld0050()  # type: ignore


def ld0074() -> tuple[list[list[float]], list[float]]:
    """生成74点 Lebedev 网格"""
    return core.ld0074()  # type: ignore


def ld0086() -> tuple[list[list[float]], list[float]]:
    """生成86点 Lebedev 网格"""
    return core.ld0086()  # type: ignore


def ld0110() -> tuple[list[list[float]], list[float]]:
    """生成110点 Lebedev 网格"""
    return core.ld0110()  # type: ignore


def ld0146() -> tuple[list[list[float]], list[float]]:
    """生成146点 Lebedev 网格"""
    return core.ld0146()  # type: ignore


def ld0170() -> tuple[list[list[float]], list[float]]:
    """生成170点 Lebedev 网格"""
    return core.ld0170()  # type: ignore


def ld0194() -> tuple[list[list[float]], list[float]]:
    """生成194点 Lebedev 网格"""
    return core.ld0194()  # type: ignore


def ld0230() -> tuple[list[list[float]], list[float]]:
    """生成230点 Lebedev 网格"""
    return core.ld0230()  # type: ignore


def ld0266() -> tuple[list[list[float]], list[float]]:
    """生成266点 Lebedev 网格"""
    return core.ld0266()  # type: ignore


def ld0302() -> tuple[list[list[float]], list[float]]:
    """生成302点 Lebedev 网格"""
    return core.ld0302()  # type: ignore


def ld0350() -> tuple[list[list[float]], list[float]]:
    """生成350点 Lebedev 网格"""
    return core.ld0350()  # type: ignore


def ld0434() -> tuple[list[list[float]], list[float]]:
    """生成434点 Lebedev 网格"""
    return core.ld0434()  # type: ignore


def ld0590() -> tuple[list[list[float]], list[float]]:
    """生成590点 Lebedev 网格"""
    return core.ld0590()  # type: ignore


def ld0770() -> tuple[list[list[float]], list[float]]:
    """生成770点 Lebedev 网格"""
    return core.ld0770()  # type: ignore


def ld0974() -> tuple[list[list[float]], list[float]]:
    """生成974点 Lebedev 网格"""
    return core.ld0974()  # type: ignore


def ld1202() -> tuple[list[list[float]], list[float]]:
    """生成1202点 Lebedev 网格"""
    return core.ld1202()  # type: ignore


def ld1454() -> tuple[list[list[float]], list[float]]:
    """生成1454点 Lebedev 网格"""
    return core.ld1454()  # type: ignore


def ld1730() -> tuple[list[list[float]], list[float]]:
    """生成1730点 Lebedev 网格"""
    return core.ld1730()  # type: ignore


def ld2030() -> tuple[list[list[float]], list[float]]:
    """生成2030点 Lebedev 网格"""
    return core.ld2030()  # type: ignore


def ld2354() -> tuple[list[list[float]], list[float]]:
    """生成2354点 Lebedev 网格"""
    return core.ld2354()  # type: ignore


def ld2702() -> tuple[list[list[float]], list[float]]:
    """生成2702点 Lebedev 网格"""
    return core.ld2702()  # type: ignore


def ld3074() -> tuple[list[list[float]], list[float]]:
    """生成3074点 Lebedev 网格"""
    return core.ld3074()  # type: ignore


def ld3470() -> tuple[list[list[float]], list[float]]:
    """生成3470点 Lebedev 网格"""
    return core.ld3470()  # type: ignore


def ld3890() -> tuple[list[list[float]], list[float]]:
    """生成3890点 Lebedev 网格"""
    return core.ld3890()  # type: ignore


def ld4334() -> tuple[list[list[float]], list[float]]:
    """生成4334点 Lebedev 网格"""
    return core.ld4334()  # type: ignore


def ld4802() -> tuple[list[list[float]], list[float]]:
    """生成4802点 Lebedev 网格"""
    return core.ld4802()  # type: ignore


def ld5294() -> tuple[list[list[float]], list[float]]:
    """生成5294点 Lebedev 网格"""
    return core.ld5294()  # type: ignore


def ld5810() -> tuple[list[list[float]], list[float]]:
    """生成5810点 Lebedev 网格"""
    return core.ld5810()  # type: ignore
