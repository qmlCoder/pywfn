import re
from pathlib import Path

path=r"D:\BaiduSyncdisk\Articles\HFV\gfile\lianben\b3lyp\gjfs\f38.gjf"
text=Path(path).read_text()

finds=re.search('(-?\d) (\d)',text)
a,b=finds.groups()
print(a,b)
