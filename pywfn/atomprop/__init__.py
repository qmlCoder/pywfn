"""
该模块用来计算原子属性
Mulliken电荷
"""
from pywfn.utils import printer

import numpy as np


from pywfn.shell import Shell

def onShell(shell:Shell):
    while True:
        printer.options('原子属性',{
            '1':'原子电荷',
            '2':'原子自旋',
            '3':'原子能量',
            '4':'原子活性',
        })
        opt=input('请输入对应选项：')
        match opt:
            case '1':
                from pywfn.atomprop import charge
                caler=charge.Calculator(shell.input.Moles()[0])
                caler.onShell(shell)
            case '2':
                from pywfn.atomprop import spin
                caler=spin.Calculator(shell.input.Moles()[0])
                caler.onShell()
            case '3':
                from pywfn.atomprop import energy
                caler=energy.Calculator(shell.input.Moles()[0])
                caler.onShell()
            case '4':
                from pywfn.atomprop import activity
                mol=shell.input.Moles(tip='输入要计算活性的分子',num=1)[0]
                caler=activity.Calculator(mol)
                caler.onShell(shell)
            case _:
                break