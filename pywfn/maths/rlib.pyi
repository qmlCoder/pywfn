f1D=list[float]
i1D=list[int]
f2D=list[list[float]]
i2D=list[list[int]]

def obt_wfns(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    coefs:f2D,
    level:int
)->tuple[list[float],list[f1D],list[f2D]]:...

def mol_rhos(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    mat_c:f2D,
    level:int
)->tuple[list[float],list[f1D],list[f2D]]:...

def march_cube(
    shape:i1D,
    grids:f2D,
    value:f1D,
    isov:float
)->tuple[f2D,i2D]:...

def a2m_weits(
    iatm:int,
    atm_grids:f2D,
    atm_weits:f1D,
    atm_pos:f2D,
    atm_rad:f1D
)->f1D:...

def lag_intpol_rs(xs:f1D,ys:f1D,ts:f1D)->f1D:...

def mat_integ_rs(
    atos:i1D,
    coes:f1D,
    alps:f1D,
    lmns:i2D,
    xyzs:f2D
)->f2D:...

def nuc_potential_rs(
    qpos:f2D,
    xyzs:f2D,
    nucs:f1D,
)->f1D:...

def ele_potential_rs(
    qpos:f2D,
    grids:f2D,
    weits:f1D,
    dens:f1D,
)->f1D:...