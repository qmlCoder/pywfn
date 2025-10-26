from pywfn import core
from numpy.typing import NDArray

def gtf_integ(alps:list[float],xyzs:list[list[float]],lmns:list[list[int]])->float:
    return core.integ.gtf_integ(alps,xyzs,lmns) # type: ignore

def local_gtf_integ(alps:list[float],xyzs:list[list[float]],syms:list[str],stm1:NDArray,stm2:NDArray)->float:
    return core.integ.local_gtf_integ(alps,xyzs,syms,stm1,stm2) # type: ignore