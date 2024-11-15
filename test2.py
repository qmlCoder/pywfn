import numpy as np
import pyvista as pv

# 假设你有一个体素数据的 NumPy 数组
# 这里我们创建一个示例数组，实际应用中你应当从文件或其他数据源加载数据
voxel_data = np.random.rand(10, 10, 10)  # 示例：随机生成的体素数据
voxel_data[voxel_data < 0.5] = 0.5  # 设定一个阈值，小于该阈值的体素被认为是空的

# 创建一个 UniformGrid 对象
grid = pv.UniformGrid()
# 设置网格的尺寸，注意这里需要加1，因为维度是从0开始计数的
grid.dimensions = np.array(voxel_data.shape) + 1
# 设置网格的原点（默认是(0,0,0)）
grid.origin = (0, 0, 0)
# 设置网格间距，默认是1
grid.spacing = (1, 1, 1)

# 将体素数据附加到网格上
# 注意，体素数据应该是一维的，所以我们需要展平它
grid.point_data["values"] = voxel_data.flatten(order="F")  # Fortran-like index order

# 创建一个 Plotter 对象
plotter = pv.Plotter()

# 添加网格到 Plotter，并指定体素数据的标量
actor = plotter.add_mesh(grid, scalars="values", opacity="linear")

# 显示可视化
plotter.show()