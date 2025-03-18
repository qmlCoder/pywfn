import sys;sys.path.append("d:/code/pywfn")

from pywfn.base import Mole
from pywfn.reader import LogReader
import matplotlib.pyplot as plt


path1=rf"D:\gfile\pywfn\5d6d\ch4_5d.out"
path2=rf"D:\gfile\pywfn\5d6d\ch4_6d.out"

reader1=LogReader(path1)
reader2=LogReader(path2)

mole1=Mole(reader1)
mole2=Mole(reader2)


# delta=mole2.coefs.CM('car')-reader2.read_CMs()[-1]
# print(delta.min(),delta.max())

# delta=mole1.coefs.SM_car-mole2.coefs.SM_car
# print(delta.min(),delta.max())

# print(mole2.coefs.SM_car)
# print(mole1.coefs.SM_car_)

delta=mole2.coefs.SM_car_-mole1.coefs.SM_car
print(delta.min(),delta.max())
# plt.matshow(delta)
# plt.show()
# delta=mole2.coefs.SM_car_-reader2.read_SM()
# print(delta.min(),delta.max())

# delta=mole2.coefs.SM_raw-reader2.read_SM()
# print(delta.min(),delta.max())

# delta=mole1.coefs.SM_raw-reader1.read_SM()
# print(delta.min(),delta.max())