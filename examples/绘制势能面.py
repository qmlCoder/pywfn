import sys
sys.path.append("D:\code\pywfn")

from pywfn.tools.getPES import Tool,Block,Blocks,Mole,Moles

tool=Tool()

molsList=[
    ('A',12),
    ('T',24),
    ('C',9)
]

molsDict={m:v for m,v in molsList}

blockData:list[list[int,str,str,str,int,int]]=[
    [0,'A','T','C',None,1], # 反应物，中间体，产物，前块，后块
    [1,'C;A','T','A',0,None],
    [2,'C','T','A',0,None],
    [3,'A','T','C',2,None],
]

blocks:list[Block]=[]
for bid,rs,ts,ps,bi,ei in blockData:
    print(bid,rs,ts,ps,bi,ei)
    rms=[Mole(float(molsDict[n])) for n in rs.split(';')]
    tms=[Mole(float(molsDict[n])) for n in ts.split(';')]
    pms=[Mole(float(molsDict[n])) for n in ps.split(';')]
    molsr,molst,molsp=Moles(rms,f'{bid}{rs}'),Moles(tms,f'{bid}{ts}'),Moles(pms,f'{bid}{ps}')
    block=Block(blocks,bid,molsr,molst,molsp,bi,ei)
    blocks.append(block)

for b,block in enumerate(blocks):
    b0=blockData[b][4] # 前一个的终点
    b1=blockData[b][5] # 后一个的起点
    block.bi0=b0
    block.bi1=b1
    print(block)

# for b,block in enumerate(blocks):
#     print('block',block)
#     print('b0',block.b0)
#     print('b1',block.b1)


tool.blocks=blocks
tool.create()
tool.write()
pass