from pywfn import config
config.IF_SHELL=True
config.IF_DEBUG=__file__=='d:\code\pywfn\main.py' # 只有在开发环境才输出log信息
if __name__=='__main__':
    from pywfn.shell import Shell
    Shell().homePage()