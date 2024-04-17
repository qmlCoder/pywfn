from pathlib import Path
import re

text=Path('./atmraddens.txt').read_text()

nums:list[str]=re.findall('\d+.\d+D[+-]\d+',text)
for i in range(0,200,10):
    num=[n.replace('D','E') for n in nums[i:i+10]]
    print(','.join(num)+',')
    

lines=text.splitlines(keepends=False)
lineNums=[]
for l,line in enumerate(lines):
    if 'npt=' not in line:continue
    num=re.search('\d+',line).group()
    num=int(num)
    lineNums.append([l,num])
print(len(lineNums))

values={}
for l,t in lineNums:
    texb='\n'.join(lines[l+1:l+t+1])
    nums=re.findall('\d\.\d+',texb)
    nums=[float(n) for n in nums]
    print(len(nums),t)
    nums+=[0]*(200-len(nums))
    values[l]=nums

f=open('extra.txt','w',encoding='utf-8')
for idx,(key,value) in enumerate(values.items()):
    # print(key,len(value))]
    lstrs=[]
    f.write(f'{idx+1}:[\n')
    for i in range(0,200,10):
        vals=value[i:i+10]
        vals=[f'{v:8.6f}' for v in vals]
        lstr=','.join(vals)+','
        f.write(f'    {lstr}\n')
    f.write('],\n')
f.close()