from pywfn.writer.xyz import XyzWriter
from pywfn.writer.gjf import GjfWriter
from pywfn.writer.cub import CubWriter
from pywfn.writer.mol import MolWriter


from pywfn.utils import printer
from pywfn.shell import Shell

def onShell(shell: Shell):
    printer.options('导出文件',{
        '1':'cub文件',
        '2':'gjf文件',
        '3':'mol文件',
        '4':'xyz文件'
    })
    opt=input('请选择文件类型:')
    mols=shell.input.Moles()
    if opt=='1':
        for mol in mols:
            CubWriter(mol).onShell()
    elif opt=='2':
        for mol in mols:
            GjfWriter(mol).onShell()
    elif opt=='3':
        for mol in mols:
            MolWriter(mol).onShell()
    elif opt=='4':
        for mol in mols:
            XyzWriter(mol).onShell()
    else:
        return