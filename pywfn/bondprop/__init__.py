from pywfn.utils import printer
from pywfn.shell import Shell
class Caler:

    def calculate(self,idx1:int,idx2:int)->float:
        return 0.0

    def print(self,resStr:str):
        printer.res(resStr)

    def resStr(self,idx1:int,idx2:int)->str:
        res=self.calculate(idx1,idx2)
        return f'{res}'

def onShell(shell:Shell):
    from pywfn.bondprop import bondOrder
    printer.info('选择要计算的键性质：')
    printer.info('1. 键级')
    mol=shell.input.Moles()[0]
    caler=bondOrder.Calculator(mol)
    caler.onShell()