import sys
sys.path.append("D:\code\pywfn")

from pywfn.tools.getPES import Tool,Block,Mole,Moles


tool=Tool()

molsList=[
    ('A',12),
    ('T',24),
    ('C',9)
]

molsDict={m:v for m,v in molsList}

blockData:list[list[str,str,str,int,int]]=[
    ['A','T','C',None,1],
    ['C;A','T','A',0,None],
    ['C','T','A',0,None],
]

blocks:list[Block]=[]
for i,(rs,ts,ps,b,e) in enumerate(blockData):
    rms=[Mole(float(molsDict[n])) for n in rs.split(';')]
    tms=[Mole(float(molsDict[n])) for n in ts.split(';')]
    pms=[Mole(float(molsDict[n])) for n in ps.split(';')]
    blocks.append(Block(i,Moles(rms,f'{i}{rs}'),Moles(tms,f'{i}{ts}'),Moles(pms,f'{i}{ps}')))

for b,block in enumerate(blocks):
    b0=blockData[b][4] # 前一个的终点
    b1=blockData[b][3] # 后一个的起点
    if b0:block.b0=blocks[b0]
    if b1:block.b1=blocks[b1]

for block in blocks:
    print(block)

tool.blocks=blocks
tool.create()
tool.write()