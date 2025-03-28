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

# sfiles=[]
# tfiles=[]

# for paths,dirnames,filenames in os.walk('pywfn'):
#     for filename in filenames:
#         path=os.path.join(paths,filename)
#         if Path(path).suffix in ['.pyc']:continue
#         if Path(path).suffix in ['.py']:
#             tfiles.append(path)
#             sfiles.append(path)
#         if Path(path).suffix in ['.dll']:
#             sfiles.append(path)
#             tfiles.append(os.path.join('libs',filename))
#         if Path(path).suffix in ['.csv','.npy']:
#             sfiles.append(path)
#             tfiles.append(os.path.join('data',filename))

# with zipfile.ZipFile('pywfn.zip','w',zipfile.ZIP_DEFLATED) as zf:
#     for sfile,tfile in zip(sfiles,tfiles):
#         zf.write(sfile,tfile)
#     zf.write('main.py','main.py')


# os.system('python -m nuitka --mingw64 --standalone --show-progress --output-dir=out --include-data-dir=pywfn/data=./data main.py')
# os.system('python -m nuitka --mingw64 --standalone --show-progress --output-dir=out main.py')
# shutil.copytree('./data','./out/main.dist/data')
# shutil.copytree('./libs','./out/main.dist/libs')

import py7zr
from datetime import datetime

params=[
    '-F',
    'main.py',
]
pyi.run(params)
if not Path('./dist/libs').exists():
    os.mkdir('./dist/libs')
shutil.copyfile('./pywfn/maths/flib.dll','./dist/libs/flib.dll')

now=datetime.now()
Y = now.year-2000
M = now.month
D = now.day
h = now.hour
m = now.minute
print(Y,M,D,h,m)
# 创建7z压缩文件
zfile=py7zr.SevenZipFile(f'pywfn_{Y}.{M:0>2}.{D:0>2}-{h:0>2}.{m:0>2}.7z', 'w')
# 遍历源文件夹中的所有文件和子目录
for root, dirs, files in os.walk('dist'):
    # 计算相对于源文件夹的相对路径
    relative_root = os.path.relpath(root, 'dist')
    for file in files:
        file_path = os.path.join(root, file)
        if '__pycache__' in file_path:continue
        if '.f90' in file_path:continue
        if '.mod' in file_path:continue
        print(file_path)
        zfile.write(file_path, os.path.join(relative_root, file))
zfile.close()