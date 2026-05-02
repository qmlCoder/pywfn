from pywfn import _core


def march_cube_v1(shape: list[int], grids, values, isov: float):
    verts, types = _core.march.march_cube_v1(shape, grids, values, isov)  # type: ignore
    return verts, types
