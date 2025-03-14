import sys;sys.path.append("d:/code/pywfn")

from pywfn.maths import flib,rlib
import numpy as np
import matplotlib.pyplot as plt

xs=np.linspace(0,4,10)
ys=np.sin(xs)
ts=np.linspace(0,4,100)

vs1=flib.lagIntpol(xs,ys,ts)
vs2=rlib.lag_intpol_rs(xs.tolist(),ys.tolist(),ts.tolist()) # type: ignore
vs2=np.array(vs2)
print(vs1)
print(vs2)

plt.scatter(xs,ys)
plt.plot(ts,vs1)
plt.plot(ts,vs2)
plt.show()