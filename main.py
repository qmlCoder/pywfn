import os
import shutil
from pathlib import Path
from pywfn import config
config.IF_SHELL=True
config.IF_DEBUG=(__file__==RF'd:\code\pywfn\main.py') # 只有在开发环境才输出log信息
# 以程序方式运行时，需要将数据文件拷贝出来，并改变默认的加载位置
try:
    shutil.copytree('pywfn/data','data')
    shutil.copytree('pywfn/libs','libs')
except:
    pass

config.ROOT_DATA=Path.cwd()/'data'
config.ROOT_LIBS=Path.cwd()/'libs'
print(f'{config.ROOT_DATA}')
print(f'{config.ROOT_LIBS}')
os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录
if __name__=='__main__':
    from pywfn.shell import Shell
    Shell().homePage()