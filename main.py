import os
import shutil
from pathlib import Path
from pywfn import config
import sys,traceback
# config.ROOT_LIBS=Path.cwd()/'libs'
os.add_dll_directory(rf"{config.ROOT_LIBS}") # 添加动态链接库目录

def error():
    exc_type,exc_value,exc_tb=sys.exc_info()
    errors=traceback.format_exception(exc_type,exc_value,exc_tb)
    error_msg = "\n".join(errors)
    print(f"程序发生错误：\n{error_msg}")
    input("按回车键退出程序...")

sys.excepthook=lambda *args:error()

if __name__=='__main__':
    from pywfn.cli import Shell
    Shell().homePage()