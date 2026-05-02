
root="e:/code/pywfn"

from pywfn.reader import CubReader
def test_cubReader():
    path=f"{root}/tests/moles/C6H6.cub"
    reader=CubReader.from_path(path)
    atoms=reader.get_atoms()
    basis=reader.get_basis()
    coefs=reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)

from pywfn.reader import GjfReader
def test_gjfReader():
    path=f"{root}/tests/moles/C2H4.gjf"
    reader=GjfReader.from_path(path)
    atoms=reader.get_atoms()
    basis=reader.get_basis()
    coefs=reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)

from pywfn.reader import LogReader
def test_logReader():
    path=f"{root}/tests/moles/C2H4.out"
    reader=LogReader.from_path(path)
    atoms=reader.get_atoms()
    basis=reader.get_basis()
    coefs=reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)

test_cubReader()
test_gjfReader()
test_logReader()