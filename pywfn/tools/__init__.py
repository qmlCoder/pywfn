"""
文件信息处理的各种工具
有需求就会有工具
"""


from pywfn import shell
from pywfn.utils import printer
def onShell(shell:"shell.Shell"):
    printer.options('实用工具',{
        '1':'分隔SCAN文件',
        '2':'分隔IRC文件',
        '3':'分隔link任务',
        '4':'提取SI信息',
    })
    opt=input('请输入选项：')
    if opt=='1':
        from pywfn.tools import spiltScan
        path=input('请输入文件路径：')
        spiltScan.Tool(path).split()
    elif opt=='2':
        from pywfn.tools import splitIrc
        path=input('请输入文件路径：')
        splitIrc.Tool().split()
    elif opt=='3':
        from pywfn.tools import splitLink
        path=input('请输入文件路径：')
        splitLink.Tool().split()
    elif opt=='4':
        from pywfn.tools import extractSI
        path=input('请输入文件路径：')
        extractSI.Tool().save()