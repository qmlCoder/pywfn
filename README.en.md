# `pywfn` -- A Python-based Wavefunction Analysis Tool

Document： https://www.xiaofei911.top/mkdocs/pywfn/


## Dependencies
```
numpy>=2.1.1
rich>=13.8.0
matplotlib>=3.9.2
```
## Usage(CLI)
``` shell
python main.py
```

## Example(API)
```python
from pywfn.base import Mol
from pywfn.reader import LogReader
from pywfn.atomProp import atomCharge

path="D:\BaiduSyncdisk\gfile\CnHn\C6H6.log" # Path to Gaussian output file
reader=LogReader(path) # Initialize log file reader 
mol=Mol(reader) # Create molecule object

caler=atomCharge.Calculator(mol) # Initialize atomic charge calculator
result=caler.mulliken() # Calculate Mulliken charges
print(result) # Print results
```

## Functions
![](./docs/pywfn_xmind.png)