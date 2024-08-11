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
    from pywfn.bondprop import order
    printer.options('键の属性',{
      '1':'各种键级'
    })
    opt=input('请输入要计算的键属性:')
    if opt=='1':
        mols=shell.input.Moles()
        caler=order.Calculator(mols[0])
        caler.onShell()