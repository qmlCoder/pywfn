import numpy as np

numDict={
    '两个苯':[0.5396]+[0.5425]*4+[0.7408]*4+[0.5650]*2,
    '三个苯':[0.6183]*4+[0.4914]*2+[0.5030]*4+[0.7702]*4+[0.5210]*2,
    '连烯烃':[0.8414]+[0.3994]*2+[0.9055]*2
}

nums=numDict['三个苯']

nums=np.array(nums)
mean=np.mean(nums)

sd=np.mean((nums-mean)**2)**0.5
print(sd)

c0=np.array([0,0,255])
c1=np.array([255,0,0])

dc=c1-c0

for num in nums:
    r,g,b=c0+num*dc
    print(f'{num:>5.2f}->{r:>3.0f},{g:>3.0f},{b:>3.0f}')
