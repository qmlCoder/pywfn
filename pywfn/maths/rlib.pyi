import numpy as np

f1D=list[float]|np.ndarray
i1D=list[int]|np.ndarray
f2D=list[list[float]]|np.ndarray
i2D=list[list[int]]|np.ndarray
f3D=list[list[list[float]]]|np.ndarray
i3D=list[list[list[int]]]|np.ndarray

def ato_wfns_rs(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    level:int
)->tuple[list[float],list[f1D],list[f2D]]:...

def obt_wfns_rs(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    coefs:f2D,
    level:int
)->tuple[list[float],list[f1D],list[f2D]]:...

def mol_rhos_rs(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    mat_c:f2D,
    level:int
)->tuple[list[float],list[f1D],list[f2D]]:...

def ato_rhos_rs(
    grids:f2D,
    xyzs:f2D,
    lmns:i2D,
    coes:f2D,
    alps:f2D,
    matp:f2D,
    level:int
)->tuple[list[f1D],list[f2D],list[f3D]]:...

def march_cube_rs(
    shape:i1D,
    grids:f2D,
    value:f1D,
    isov:float
)->tuple[f2D,i2D,f2D,i2D]:...

def a2m_weits_rs(
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

def get_grids_rs(nx:int,ny:int,nz:int)->f2D:...

def ele_mat_rs(matc:f2D,mats:f2D)->f2D:...

def calc_EIEBA_rs(
    fragA:f2D,
    fragB:f2D,
    rhosA:f1D,
    rhosB:f1D,
    gridsA:f2D,
    weitsA:f1D,
    gridsB:f2D,
    weitsB:f1D,
)->tuple[float,float,float]:...

def search_title(path:str)->dict[str,int]:...

def read_geome(path:str,title_number:int)->tuple[f1D,f2D]:...
def read_nbas(path:str,title_number:int)->int:...

def read_cmat(path:str,title_number:int,nbas:int):...
def read_basis(path:str,title_number:int)->list[tuple[int,int,int,float,float]]:...