import os
import shutil
from pathlib import Path
from pywfn import config
config.IF_SHELL=True
config.IF_DEBUG=(__file__==RF'd:\code\pywfn\main.py') # 只有在开发环境才输出log信息
# 以程序方式运行时，需要将数据文件拷贝出来，并改变默认的加载位置
if Path('pywfn').exists():
    # if not Path('data').exists():Path('data').mkdir()
    if not Path('libs').exists():Path('libs').mkdir()
    cwd=Path().cwd()
    for root,dirs,files in os.walk(f'{cwd}/pywfn/data'):
        for name in files:
            if Path(name).suffix in ['.pyc','.py']:continue
            shutil.copyfile(f'{root}/{name}',f'data/{name}')
            
    for root,dirs,files in os.walk(f'{cwd}/pywfn/libs'):
        for name in files:
            if Path(name).suffix in ['.f90','.mod','.o']:continue
            shutil.copyfile(f'{root}/{name}',f'libs/{name}')


config.ROOT_DATA=Path.cwd()/'data'
config.ROOT_LIBS=Path.cwd()/'libs'
# print(f'{config.ROOT_DATA}')
# print(f'{config.ROOT_LIBS}')
# os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录
if __name__=='__main__':
    from pywfn.shell import Shell
    Shell().homePage()