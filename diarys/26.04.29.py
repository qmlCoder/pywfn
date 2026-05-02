"""
mocv法指定原子轨道绘制电子密度
"""
from pywfn.base import Mole
from pywfn.orbtprop import obtmat
from pywfn.gridprop import density,CubeGrid
from pywfn.atomprop import direction
from pywfn.orbtprop.obtmat import Mocv
from pywfn.writer import CubWriter

# 加载分子
mole = Mole.from_file(rf"c:\Users\11032\Desktop\gfile\pywfn\C6H6.fch")

dir_caler=direction.Calculator(mole) #方向计算器
natm=mole.atoms.len() # 原子个数
mocvs={} # mocv方法用到的参数
for atm in range(natm):
    stm=dir_caler.LCS(atm) # 原子的局部坐标系
    mocv=Mocv.new(stm,[0],[1,1,0],[0,0,0,0,0]) # s,p,d原子轨道保留的成分
    mocvs[atm]=mocv
cmat_caler = obtmat.Calculator(mole) #分子轨道计算器
cmat_mocv=cmat_caler.mocv(mocvs) #计算mocv方法的分子轨道
mole.set_cmat("sph",cmat_mocv) #修改分子轨道为mocv方法的分子轨道

dens_caler=density.Calculator(mole) #电子密度计算器
p0,p1=mole.border() # 分子边界（包围分子的盒子）
grid_caler=CubeGrid() #格点计算器
step=0.2 #格点步长
bord=4.0 #格点扩展
grid_caler.set_v1(p0,p1,step,bord) # 设置格点参数
shape,grids=grid_caler.get() # 格点形状和数值
dens=dens_caler.mol_rho_cm(grids,0)[0].reshape(1,-1) # 计算电子密度 [格点数,轨道数]

writer=CubWriter() # cub文件写入器
writer.read_mole(mole) # 从分子中读取原子类型和坐标
x0,y0,z0=p0
writer.set(
    obts=[-1], # 因为是电子密度，所以轨道序号为-1
    pos0=[x0-bord,y0-bord,z0-bord], # 真实格点盒子比分子盒子大
    size=shape, # 格点形状
    step=[step,step,step], # 格点步长
    vals=dens # 格点数值（电子密度）
)
writer.save(rf"E:\code\pywfn\tests\moles\sigma_rho.cub") # 保存为cub文件