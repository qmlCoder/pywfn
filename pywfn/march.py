from pywfn import core

def march_cube_v1(shape:list[int],grids,values,isov:float):
    verts,types=core.march.march_cube_v1(shape,grids,values,isov) # type: ignore
    return verts,types
