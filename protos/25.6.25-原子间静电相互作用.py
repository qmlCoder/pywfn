import sys;sys.path.append('d:/code/pywfn')

from pywfn.base import Mole
from pywfn.reader import LogReader
from pywfn.gridprop import density,dftgrid
from pywfn.fragprop import energy
from pywfn.maths import rlib

import numpy as np
import matplotlib.pyplot as plt

# 设置中文字体和解决负号显示问题
plt.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体或其他中文字体
plt.rcParams['axes.unicode_minus'] = False    # 解决负号 '-' 显示为方块的问题

files=['08','09','10','11','12','13',]
eng_ees=[]
eng_nes=[]
eng_nns=[]
mol_engs=[]
for file in files:
    mole=Mole(LogReader(rf"d:\gfile\calc\S{file}.out"))
    caler=energy.Calculator(mole)
    eng_ee,eng_ne,eng_nn=caler.EIEBA([1],[2])
    eng_ees.append(eng_ee)
    eng_nes.append(eng_ne)
    eng_nns.append(eng_nn)
    mol_engs.append(mole.energy)

eng_ees=np.array(eng_ees)
eng_nes=np.array(eng_nes)
eng_nns=np.array(eng_nns)

fig,axs=plt.subplots(1,5)
xs=np.linspace(0.8,1.3,6)
axs[0].plot(xs,eng_ees)
axs[1].plot(xs,eng_nes)
axs[2].plot(xs,eng_nns)
axs[3].plot(xs,eng_ees+eng_nes+eng_nns)
axs[4].plot(xs,mol_engs)

axs[0].set_title('电子-电子排斥')
axs[1].set_title('核-电子吸引')
axs[2].set_title('核-核排斥')
axs[3].set_title('静电相互作用')
axs[4].set_title('分子能量')
plt.show()