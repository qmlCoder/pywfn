from setuptools import setup,find_packages
from pathlib import Path
from pywfn import __version__

HERE=Path(__file__).parent
README=(HERE / "README.md").read_text()

setup(
    name="pywfn",
    version=__version__,
    description="A Python library for wave functions analysis",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://www.xiaofei911.top/mkdocs/pywfn/",
    author='Xiaofei Shi',
    packages=find_packages(),
    include_package_data=True,
    install_requires=[],
    entry_points={
        
    }
)