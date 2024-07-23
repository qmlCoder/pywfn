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

BASE_VECTOR=np.array([0,0,1]) # 标准向量，求原子法向量和轨道方向的时候要与该向量夹角小于90°

IF_ORBITAL_ORDER=False #是否计算键级中每个轨道的成分，尤其是mayer键级拆分成每个轨道的成分很麻烦

RENDER_ATOM_RANGE=20000  # 原子轨道的渲染范围

RENDER_CLOUD_STEP=get_config('render.cloud.step',0.2)  # 导出cub文件时的步长

RENDER_CLOUD_BORDER=get_config('render.cloud.border',4.0) # 到处cub文件的边框

IMG_SPACE_STEP=get_config('img.space.step',0.1)

BOHR_RADIUS=1/0.529177

SHOW_PRINT=False

IF_CM_P=False # 是否使用投影轨道

ROOT_DATA=Path(__file__).parent/'data'
ROOT_LIBS=Path(__file__).parent/'libs'
# os.add_dll_directory(rf"{ROOT_LIBS}") # 添加动态链接库目录