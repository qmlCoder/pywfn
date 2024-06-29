"""
文件信息处理的各种工具
有需求就会有工具

工具应该只是尽可能只提取信息，或需能提取信息
信息到文件的步骤在shell层完成
"""

from pywfn import shell
from pywfn.utils import printer
def onShell(shell:"shell.Shell"):
    printer.options('实用工具',{
        '1':'分隔SCAN文件',
        '2':'分隔IRC 文件',
        '3':'分隔link任务',
        '4':'提取 SI 信息',
    })
    opt=input('请输入选项：')
    if opt=='1': # 分隔SCAN文件
        from pywfn.tools import spiltScan
        path=input('请输入文件路径：')
        spiltScan.Tool(path).save()
    elif opt=='2': # 分隔IRC文件
        from pywfn.tools import splitIrc
        path=input('请输入文件路径：')
        splitIrc.Tool(path).split()
    elif opt=='3': # 分隔link任务
        from pywfn.tools import splitLink
        paths=shell.input.Moles(mtype='path')
        for path in paths:
            splitLink.Tool(path).split()
    elif opt=='4': # 提取SI信息
        from pywfn.tools import extractSI
        from pywfn.utils import parse_intList
        print('(1:能量,2:坐标,3:频率):')
        msgStr=input('请输入需要保存的信息[*]')
        if msgStr:
            msgs=parse_intList(msgStr)
        else:
            msgs=[1,2,3]
        same=shell.input.Bool('是否保存到同一文件内？',default=True)
        mols=shell.input.Moles()
        for mol in mols:
            tool=extractSI.Tool(mol)
            tool.selects=msgs
            tool.sameFile=same
            tool.save()