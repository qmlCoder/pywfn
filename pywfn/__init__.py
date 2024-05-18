"""
基于python的波函数分析程序
"""
import numpy as np
np.set_printoptions(precision=4, suppress=True)

# __version__='0.0.1'
# __description__='A wave function analysis program based on python'
# __author__='XiaoFei Shi'
# __email__='1103275712@qq.com'
# __url__='https://www.xiaofei911.top/mkdocs/pywfn/'

def main():
    from pywfn.shell import Shell
    Shell().homePage()