"""
该脚本将一些非py文件拷贝到根目录
"""
from pathlib import Path
import zipfile
import os
import shutil
import os
import shutil
from PyInstaller import __main__ as pyi

# os.system('python -m nuitka --mingw64 --standalone --show-progress --output-dir=out  --include-data-dir=pywfn/maths=./libs main.py')


import py7zr
from datetime import datetime

# params=[
#     '-F',
#     'main.py',
# ]
# pyi.run(params)

now=datetime.now()
Y = now.year-2000
M = now.month
D = now.day
h = now.hour
m = now.minute
# print(Y,M,D,h,m)
# 创建7z压缩文件
zfile=py7zr.SevenZipFile(f'pywfn_{Y}.{M:0>2}.{D:0>2}-{h:0>2}.{m:0>2}.7z', 'w')
# 遍历源文件夹中的所有文件和子目录
for root, dirs, files in os.walk(rf'out\main.dist'):
    # 计算相对于源文件夹的相对路径
    relative_root = os.path.relpath(root, rf'out\main.dist')
    for file in files:
        file_path = os.path.join(root, file)
        if '__pycache__' in file_path:continue
        if '.f90' in file_path:continue
        if '.mod' in file_path:continue
        print(file_path)
        zfile.write(file_path, os.path.join(relative_root, file))
zfile.close()