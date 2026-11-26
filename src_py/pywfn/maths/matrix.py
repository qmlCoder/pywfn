from pywfn import _core


def ele_mat(matc: list[list[float]], mats: list[list[float]]) -> list[list[float]]:
    """计算矩阵的逐元素乘积（Hadamard 积）

    Args:
        matc: 第一个矩阵
        mats: 第二个矩阵

    Returns:
        两个矩阵逐元素相乘的结果
    """
    return _core.maths.matrix.ele_mat(matc, mats)  # type: ignore
