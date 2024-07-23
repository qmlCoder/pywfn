"""
文件信息处理的各种工具
有需求就会有工具

工具应该只是尽可能只提取信息，或需能提取信息
信息到文件的步骤在shell层完成
"""

from pywfn import shell
from pywfn.utils import printer

from pathlib import Path
cwd=Path.cwd()

def onShell(shell:"shell.Shell"):
    printer.options('实用工具',{
        '1':'分割SCAN 文件',
        '2':'分割 IRC 文件',
        '3':'分割link 任务',
        '4':'拼接 gjf 文件',
        '5':'提取  SI 信息',
        '6':'环心添加Bq原子'
        
    })
    opt=input('请输入选项：')
    if opt=='1': # 分割SCAN文件
        from pywfn.tools import spiltScan
        paths=shell.input.Paths()
        for path in paths:
            spiltScan.Tool(path).save()
    elif opt=='2': # 分割IRC文件
        from pywfn.tools import splitIrc
        paths=shell.input.Paths()
        for path in paths:
            splitIrc.Tool(path).split()
    elif opt=='3': # 分割link任务
        from pywfn.tools import splitLink
        paths=shell.input.Paths()
        for path in paths:
            splitLink.Tool(path).split()
    elif opt=='4':
        from pywfn.tools import joinGjf
        paths=shell.input.Paths()
        joinGjf.Tool(paths).save(f'{cwd}/join.gjf')
    elif opt=='5': # 提取SI信息
        from pywfn.tools import extractSI
        paths=shell.input.Paths()
        tool=extractSI.Tool(paths)
        path0=Path(paths[0]) # 第一个分子
        spath=path0.parent/'SI.txt' # 第一个文件所在的文件夹
        tool.save(f'{spath}')
    elif opt=='6':
        from pywfn.tools import ringBq
        printer.info('在gjf文件指定环的中心添加Bq原子，方便NICS计算')
        mol=shell.input.Moles(count=1)[0]
        tool=ringBq.Tool(mol)
        tool.onShell(shell)