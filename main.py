from pywfn import config
from pywfn.shell import Shell
if __name__=='__main__':
    config.IF_SHELL=True
    Shell().homePage()