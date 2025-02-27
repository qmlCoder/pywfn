# 用来管理程序的一些设置
from pathlib import Path
import numpy as np
import json
import os

configPath=Path.cwd()/'config.json'
if not configPath.exists():
    configPath.touch()
    configPath.write_text('{}')

configDict:dict=json.loads(configPath.read_text())

def get_config(key:str,value):
    if key not in configDict.keys():
        configDict[key]=value
        set_config(key,value)
    return configDict[key]

def set_config(key:str,value):
    global configDict
    configDict[key]=value
    configPath.write_text(json.dumps(configDict,indent=1))

BOND_LIMIT=get_config('bond.limit',1.8*1.889) # 判断量原子之间是否成键的长度限制
IF_DEBUG=True # 是否开启debug,控制项目中所有的打印,避免与shell中的print冲突
IF_SHELL=False # 是否在shell中执行
IF_BUFF=get_config('if_buff',True) # 是否生成缓存

# BASE_VECTOR=np.array([1,0,0]) # 标准向量，求原子法向量和轨道方向的时候要与该向量夹角小于90°
BASE_VECTOR=np.random.rand(3) # 标准向量，求原子法向量和轨道方向的时候要与该向量夹角小于90°
BASE_VECTOR/=np.linalg.norm(BASE_VECTOR)


IF_ORBITAL_ORDER=False #是否计算键级中每个轨道的成分，尤其是mayer键级拆分成每个轨道的成分很麻烦

RENDER_ATOM_RANGE=20000  # 原子轨道的渲染范围

RENDER_CLOUD_STEP=get_config('render.cloud.step',0.2)  # 导出cub文件时的步长

RENDER_CLOUD_BORDER=get_config('render.cloud.border',4.0) # 到处cub文件的边框

IMG_SPACE_STEP=get_config('img.space.step',0.1)

BOHR_RADIUS=1/0.529177

SHOW_PRINT=False

IF_CM_P=False # 是否使用投影轨道
IF_BUFF=get_config('if_buff',True) # 是否生成缓存

# ROOT_DATA=Path(__file__).parent/'data'
ROOT_LIBS=Path(__file__).parent/'maths'
GJF_TITLE=get_config('gjf.title','b3lyp/6-31g(d) pop=full gfinput iop(3/33=1)')

MARCH_ISOV_WFNS = get_config('march.isov',0.03) # Marching Cubes算法波函数等值面阈值
MARCH_ISOV_DENS = get_config('march.dens',0.001) # Marching Cubes算法电子密度等值面阈值

# os.add_dll_directory(rf"{ROOT_LIBS}") # 添加动态链接库目录

CACHE=False # 是否使用缓存