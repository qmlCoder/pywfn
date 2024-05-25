import sys
sys.path.append("D:\code\pywfn")

from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.bondprop import bondOrder

import numpy as np
import matplotlib.pyplot as plt

root='D:/BaiduSyncdisk/gfile/scans/lianxi/lianxiScan'

orderList=[]
angleList=[]
energys=[]
for i in range(37):
    path=f'{root}/f{i+1:0>2}.log'
    mol=Mol(LogReader(path))
    caler=bondOrder.Calculator(mol)
    # orders=caler.dirMayer([[1,5],[1,3],[3,8]])
    # print(orders.shape)
    # if orders.shape[0]==6:
    #     order=orders[[1,3,5],-1]
    # else:
    #     order=orders[:,-1]
    orders=caler.piOrder()
    order=orders[[1,2,4],-1]
    angle=mol.params([5,1,3,8])
    angleList.append(angle)
    orderList.append(order)
    energys.append(mol.energy)
orderList=np.array(orderList)
print(orderList.shape)
means=np.mean(orderList,axis=1).reshape(-1,1)
print(means.shape)
stds=np.std(orderList,axis=1)
print(stds.shape)
angs=np.linspace(-np.pi,np.pi,37)

fig,axs=plt.subplots(1,4,figsize=(16,4))
axs[0].plot(angs,orderList[:,0],marker='.')
axs[0].plot(angs,orderList[:,1],marker='.')
axs[0].plot(angs,orderList[:,2],marker='.')

means=means.flatten()
for i in np.linspace(0,1,10):
    vals=-(1-i)*stds+i*means # 标准差越小，芳香性越强，方差越大，芳香性越强
    
    axs[1].plot(angs,vals,marker='.')
    vals=(means**i)/(stds**(1-i))
    axs[2].plot(angs,vals,marker='.')
axs[3].plot(angs,-np.array(energys),'.-b')
plt.show()