"""
该脚本将一些非py文件拷贝到根目录
"""
from pathlib import Path
import zipfile
import os

sfiles=[]
tfiles=[]

for paths,dirnames,filenames in os.walk('pywfn'):
    for filename in filenames:
        path=os.path.join(paths,filename)
        if Path(path).suffix in ['.pyc']:continue
        if Path(path).suffix in ['.py']:
            tfiles.append(path)
            sfiles.append(path)
        if Path(path).suffix in ['.dll']:
            sfiles.append(path)
            tfiles.append(os.path.join('dlls',filename))
        if Path(path).suffix in ['.csv','.npy']:
            sfiles.append(path)
            tfiles.append(os.path.join('data',filename))

with zipfile.ZipFile('pywfn.zip','w',zipfile.ZIP_DEFLATED) as zf:
    for sfile,tfile in zip(sfiles,tfiles):
        zf.write(sfile,tfile)
    zf.write('main.py','main.py')