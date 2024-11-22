from pywfn.writer.xyz import XyzWriter
from pywfn.writer.gjf import GjfWriter
from pywfn.writer.cub import CubWriter
from pywfn.writer.mol import MolWriter
from pywfn.writer.obj import ObjWriter
from pywfn.utils import printer
from pywfn.shell import Shell

from pathlib import Path
from pywfn.data import temps

def onShell(shell: Shell):
    printer.options('文件转换',{
        '0':'导出模板',
        '1':'gjf文件',
        '2':'mol文件',
        '3':'xyz文件',
    })
    opt=input('请选择文件类型:')
    
    match opt:
        case '0':
            (Path().cwd()/'temp.gjf.txt').write_text(temps.gjf)
            (Path().cwd()/'temp.mol.txt').write_text(temps.mol)
            (Path().cwd()/'temp.mol.xyz').write_text(temps.xyz)
        case '1':
            mols=shell.input.Moles()
            for mol in mols:
                writer=GjfWriter().fromMol(mol)
                path=Path(mol.reader.path).with_suffix('.gjf')
                writer.save(f'{path}')
        case '2':
            mols=shell.input.Moles()
            for mol in mols:
                writer=MolWriter().fromMol(mol)
                path=Path(mol.reader.path).with_suffix('.mol')
                writer.save(f'{path}')
        case '3':
            mols=shell.input.Moles()
            for mol in mols:
                writer=XyzWriter().fromMol(mol)
                path=Path(mol.reader.path).with_suffix('.xyz')
                writer.save(f'{path}')
        case _:
            return