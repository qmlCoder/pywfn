"""
各种实用工具的合集
调用pywfn核心函数，或者是对文件进行处理
- editGjf   编辑gjf文件
- engCorr   能量矫正
- extractSI 提取SI信息
- getPES    绘制势能面
- joinGjf   拼接gjf文件
- ringBq    环心添加Bq原子
- scanGjf   生成刚性扫描的一堆gjf文件
- splitIrc  分割IRC 文件
- splitLink 分割link任务
- spiltScan 分割SCAN文件
- xmlEdit   编辑xml文件
- xmlShow   显示xml文件为图像
"""

from pywfn import shell
from pywfn.utils import printer

from pathlib import Path
cwd=Path.cwd()

def onShell(shell:"shell.Shell"):
    printer.options('实用工具',{
        '1':'提取  SI 信息',
        '2':'分割SCAN 文件',
        '3':'分割 IRC 文件',
        '4':'分割link 任务',
        '5':'拼接 gjf 文件',
        '6':'环心添加Bq原子',
        '7':'加/减gjf电子数'
        
    })
    opt=input('请输入选项：')
    match opt:
        case '1': # 提取SI信息
            from pywfn.tools import extractSI
            paths=shell.input.Paths()
            tool=extractSI.Tool(paths)
            path0=Path(paths[0]) # 第一个分子
            spath=path0.parent/'SI.txt' # 第一个文件所在的文件夹
            tool.save(f'{spath}')
        case '2': # 分割SCAN文件
            from pywfn.tools import spiltScan
            paths=shell.input.Paths()
            for path in paths:
                spiltScan.Tool(path).save()
        case '3': # 分割IRC文件
            from pywfn.tools import splitIrc
            paths=shell.input.Paths()
            for path in paths:
                splitIrc.Tool(path).split()
        case '4': # 分割link任务
            from pywfn.tools import splitLink
            paths=shell.input.Paths()
            for path in paths:
                splitLink.Tool(path).split()
        case '5': # 拼接gjf文件
            from pywfn.tools import joinGjf
            paths=shell.input.Paths()
            joinGjf.Tool(paths).save(f'{cwd}/join.gjf')
        case '6': # 环心添加Bq原子
            from pywfn.tools import ringBq
            printer.info('在gjf文件指定环的中心添加Bq原子，方便NICS计算')
            mol=shell.input.Moles(count=1)[0]
            tool=ringBq.Tool(mol)
            tool.onShell(shell)
        case '7': # 加/减gjf电子数
            from pywfn.tools import editGjf
            mols=shell.input.Moles()
            nele=input('输入加减电子数: ')
            nele=int(nele)
            if nele>0:
                sufx='n'
            elif nele<0:
                sufx='p'
            else:
                raise ValueError('加减电子数不能为0')
            for mol in mols:
                tool=editGjf.Tool(mol)
                oldCharge,oldSpin=tool.addEle(nele)
                path=Path(mol.reader.path)
                tool.save(f'{path.parent}/{path.stem}_{sufx}.gjf')
                mol.props.set('charge',oldCharge)
                mol.props.set('spin',oldSpin)