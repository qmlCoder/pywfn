from pathlib import Path

from pywfn.reader import CubReader, GjfReader, LogReader, get_reader

root = Path(__file__).parent.parent


def test_cubReader():
    path = f"{root}/moles/C6H6.cub"
    reader = CubReader.from_path(path)
    atoms = reader.get_atoms()
    basis = reader.get_basis()
    coefs = reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)


def test_gjfReader():
    path = f"{root}/moles/C2H4.gjf"
    reader = GjfReader.from_path(path)
    atoms = reader.get_atoms()
    basis = reader.get_basis()
    coefs = reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)


def test_logReader():
    path = f"{root}/moles/C2H4.out"
    reader = LogReader.from_path(path)
    atoms = reader.get_atoms()
    basis = reader.get_basis()
    coefs = reader.get_coefs()
    print(atoms)
    print(basis)
    print(coefs)


def test_get_reader():
    """按后缀自动选择读取器，并保证返回的实例可用"""
    cases = {
        "C6H6.fch": "FchReader",
        "C2H4.out": "LogReader",
        "C2H4.gjf": "GjfReader",
        "C2H4.cub": "CubReader",
    }
    for name, expect in cases.items():
        reader = get_reader(f"{root}/moles/{name}")
        assert type(reader).__name__ == expect, name
        # reader 返回的是 core 层的 Atoms，`syms` 在那里是方法
        atoms = reader.get_atoms()
        assert atoms is not None and not isinstance(atoms, str), name
        syms = atoms.syms() if callable(atoms.syms) else atoms.syms
        assert len(syms) > 0, name

    # 不支持的后缀应当抛出 ValueError，而不是把路径当成 core 对象
    try:
        get_reader(f"{root}/moles/not_exists.abc")
    except ValueError:
        pass
    else:
        raise AssertionError("不支持的后缀应当抛出 ValueError")


if __name__ == "__main__":
    test_cubReader()
    test_gjfReader()
    test_logReader()
    test_get_reader()
