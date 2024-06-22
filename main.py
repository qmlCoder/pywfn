import os
# os.add_dll_directory(rf"D:\program\mingw64\bin")
os.add_dll_directory(rf"D:\code\pywfn\libs")

from pathlib import Path
from pywfn import config
config.IF_SHELL=True
config.IF_DEBUG=(__file__==RF'd:\code\pywfn\main.py') # 只有在开发环境才输出log信息
# 以程序方式运行时，需要将数据文件拷贝出来，并改变默认的加载位置
config.ROOT_DATA=Path.cwd()/'data'
config.ROOT_LIBS=Path.cwd()/'libs'
print(f'{config.ROOT_DATA}')
print(f'{config.ROOT_LIBS}')

if __name__=='__main__':
    from pywfn.shell import Shell
    Shell().homePage()