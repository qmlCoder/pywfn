from pywfn.writer.xyz import XyzWriter
from pywfn.writer.gjf import GjfWriter
from pywfn.writer.cub import CubWriter
from pywfn.writer.mol import MolWriter
from pywfn.writer.obj import ObjWriter


from pywfn.utils import printer
from pywfn.shell import Shell

def onShell(shell: Shell):
    printer.options('文件转换',{
        '1':'gjf文件',
        '2':'mol文件',
        '3':'xyz文件'
    })
    opt=input('请选择文件类型:')
    mols=shell.input.Moles()
    match opt:
        case '1':
            for mol in mols:
                GjfWriter(mol).onShell()
        case '2':
            for mol in mols:
                MolWriter(mol).onShell()
        case '3':
            for mol in mols:
                XyzWriter(mol).onShell()
        case _:
            return