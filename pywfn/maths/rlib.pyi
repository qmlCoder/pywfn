f1D=list[float]
i1D=list[int]
f2D=list[list[float]]
i2D=list[list[int]]

def mol_rhos(
        grids:f2D,
        xyzs:f2D,
        lmns:i2D,
        coes:f2D,
        alps:f2D,
        mat_c:f2D,
        level:int
    )->tuple[float,f1D,f2D]:...

def march_cube(
    shape:i1D,
    grids:f2D,
    value:f1D,
    isov:float
)->tuple[tuple[float]]:...