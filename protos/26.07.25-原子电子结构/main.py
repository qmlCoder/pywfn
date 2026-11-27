"""
绘制原子电子结构的图像
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pywfn.datas import elems



boxs=[
    [1,'s',1,0],
    [2,'s',1,24.5],
    [3,'s',1,49],
    [4,'s',1,73.5],
    [5,'s',1,98],
    [6,'s',1,122.5],
    [7,'s',1,147],
    [2,'p',3,40],
    [3,'p',3,64.5],
    [4,'p',3,89],
    [5,'p',3,113.5],
    [6,'p',3,138],
    [7,'p',3,162.5],
    [3,'d',7,84],
    [4,'d',7,108.5],
    [5,'d',7,133],
    [6,'d',7,157.5],
    [4,'f',13,128],
    [5,'f',13,152.5],
]
txts=[
    [1,147,0],
    [3,162.5,-1],
    [4,162.5,0],
    [5,162.5,1],
    [7,157.5,-2],
    [8,157.5,-1],
    [9,157.5,0],
    [10,157.5,1],
    [11,157.5,2],
    [13,152.5,-3],
    [14,152.5,-2],
    [15,152.5,-1],
    [16,152.5,0],
    [17,152.5,1],
    [18,152.5,2],
    [19,152.5,3],

]
boxs.sort(key=lambda x: x[-1])
nums={'s': 1, 'p': 3, 'd': 5, 'f': 7}

lasts=[0,2,10,18,36,54,86,118]
# 获取每个元素的内层轨道数
def get_inner(idx):
    for i in range(len(lasts)-1):
        if lasts[i]<idx<=lasts[i+1]:
            ele=lasts[i]
            return ele//2 + (0 if ele%2==0 else 1)
    return 0

def draw_elem(idx):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 200)
    ax.set_ylim(-10, 200)
    ax.axis('off')
    inner=get_inner(idx+1) # 内层轨道数
    print(idx+1,inner)
    recs=[]
    for box in boxs:
        n, sym, x, y = box
        for j in range(nums[sym]):
            px=(x+j)*10
            recs.append([px,y])
            facecolor=(0,0,0,0.1) if len(recs)<=inner else 'none'
            rect = Rectangle((px, y), 8, 20, facecolor=facecolor, edgecolor='black')

            ax.add_patch(rect)
        ax.text(x*10-5, y + 10, f'{n}{sym}', ha='center', va='center')
    for x,y,txt in txts:
        ax.text(x*10+4, y+22, txt, ha='center', va='bottom')
    eles=elems.ELES()[idx]

    for i in range(len(eles)):
        ele=eles[i]
        px, y = recs[i]
        head=1.4
        if ele>=1:
            ax.arrow(px+2.5, y+2, 0, 14, head_width=head, head_length=head, color='black')
        if ele>=2:
            ax.arrow(px+2.5+3, y+18, 0, -14, head_width=head, head_length=head, color='black')
    fig.savefig(f'./imgs/elems_{idx+1:0>3}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

for i in range(118):
    draw_elem(i)
