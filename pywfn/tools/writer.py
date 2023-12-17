"""
需要导出文件：gjf等
"""
from pathlib import Path
from pywfn.data import temps

def gjf(path,chk:str,coord:list[list[float]],charge:int,spin:int):
    temp=temps.gjf
    replaces=[
        ['<COORD>', cord2str(coord)],
        ['<CHK>', chk],
        ['<CHARGE>',f'{charge}'],
        ['<MULTI>',f'{spin}']
    ]
    for k,v in replaces:
        temp=temp.replace(k,v)
    Path(path).write_text(temp)

def cord2str(coord:list[list[str,float,float,float]]):
    strs=[]
    for s,x,y,z in coord:
        strs.append(f' {s:<6}{x:>14.8f}{y:>14.8f}{z:>14.8f}')
    return '\n'.join(strs)