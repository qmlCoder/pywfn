"""
该模块用来计算原子属性
Mulliken电荷
"""
from pywfn.utils import printer

import numpy as np

class AtomCaler:
    def __init__(self) -> None:
        self.logTip:str=''

    def resStr(self)->np.ndarray:pass

    def logRes(self):
        printer.info(self.logTip)
        printer.res(self.resStr())


from pywfn.shell import Shell

def onShell(shell:Shell):
    printer.info('1. 原子电荷')
    printer.info('2. 原子自旋')
    printer.info('3. 原子能量')
    printer.info('4. 原子活性')
    printer.info('5. 原子活性')
    opt=printer.info('请输入对应选项：')
    if opt=='':return
    if opt=='1':
        from pywfn.atomprop import charge
        caler=charge.Calculator(shell.input.Moles()[0])
        caler.onShell()
    if opt=='2':
        from pywfn.atomprop import spin
        caler=spin.Calculator(shell.input.Moles()[0])
        caler.onShell()
    if opt=='3':
        from pywfn.atomprop import energy
        caler=energy.Calculator(shell.input.Moles()[0])
        caler.onShell()
    if opt=='4':
        from pywfn.atomprop import activity
        caler=activity.Calculator(shell.input.Moles()[0])
        caler.onShell()
