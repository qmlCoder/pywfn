"""
将python的lebedev转为Fortran的
"""

from pywfn.data import lebedev
import re
from pathlib import Path
import numpy as np

text0="""
module Lebdev
    implicit none
    integer,parameter::lens(6)=[6,12,8,24,24,48] !每个code对应的点数
contains
subroutine gen_oh(code,a,b,v,xs,ys,zs,ws,l)
    integer, intent(in) :: code
    integer, intent(in) :: l
    real(8), intent(inout) :: a,b
    real(8), intent(in) :: v
    real(8), intent(out):: xs(l),ys(l),zs(l),ws(l)
    real(8)::z=0.0d0,c
    select case(code)
    case(1)
        a=1.0d0
        xs=[ a,-a, z, z, z, z]
        ys=[ z, z, a,-a, z, z]
        zs=[ z, z, z, z, a,-a]
        ws=[ v, v, v, v, v, v]
    case(2)
        a=sqrt(1.0/2.0)
        xs=[ z, z, z, z, a,-a, a,-a, a,-a, a,-a]
        ys=[ a,-a, a,-a, z, z, z, z, a, a,-a,-a]
        zs=[ a, a,-a,-a, a, a,-a,-a, z, z, z, z]
        ws=[ v, v, v, v, v, v, v, v, v, v, v, v]
    case(3)
        a=sqrt(1.0/3.0)
        xs=[ a,-a, a,-a, a,-a, a,-a]
        ys=[ a, a,-a,-a, a, a,-a,-a]
        zs=[ a, a, a, a,-a,-a,-a,-a]
        ws=[ v, v, v, v, v, v, v, v]
    case(4)
        b = sqrt(1.0-2*a**2)
        xs=[ a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, b,-b, b,-b, b,-b, b,-b]
        ys=[ a, a,-a,-a, a, a,-a,-a, b, b,-b,-b, b, b,-b,-b, a, a,-a,-a, a, a,-a,-a]
        zs=[ b, b, b, b,-b,-b,-b,-b, a, a, a, a,-a,-a,-a,-a, a, a, a, a,-a,-a,-a,-a]
        ws=[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    case(5)
        b=sqrt(1.0-a**2)
        xs=[ a,-a, a,-a, b,-b, b,-b, a,-a, a,-a, b,-b, b,-b, z, z, z, z, z, z, z, z]
        ys=[ b, b,-b,-b, a, a,-a,-a, z, z, z, z, z, z, z, z, a,-a, a,-a, b,-b, b,-b]
        zs=[ z, z, z, z, z, z, z, z, b, b,-b,-b, a, a,-a,-a, b, b,-b,-b, a, a,-a,-a]
        ws=[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    case (6)
        c=sqrt(1.0-a**2-b**2)
        xs=[ a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, b,-b, b,-b, b,-b, b,-b, b,-b, b,-b, b,-b, b,-b, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c]
        ys=[ b, b,-b,-b, b, b,-b,-b, c, c,-c,-c, c, c,-c,-c, a, a,-a,-a, a, a,-a,-a, c, c,-c,-c, c, c,-c,-c, a, a,-a,-a, a, a,-a,-a, b, b,-b,-b, b, b,-b,-b]
        zs=[ c, c, c, c,-c,-c,-c,-c, b, b, b, b,-b,-b,-b,-b, c, c, c, c,-c,-c,-c,-c, a, a, a, a,-a,-a,-a,-a, b, b, b, b,-b,-b,-b,-b, a, a, a, a,-a,-a,-a,-a]
        ws=[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    end select
end subroutine gen_oh

subroutine ldfunc(as,bs,vs,cs,li,xs,ys,zs,ws,lo)
    integer,intent(in)::li,lo
    real(8), intent(in) :: as(li),bs(li),vs(li)
    integer,intent(in)::cs(li) !code
    real(8), intent(out) ::  xs(lo),ys(lo),zs(lo),ws(lo)
    real(8)::a,b,v
    integer::i
    integer::u,l !索引的上下界
    integer::code
    
    integer::idx=1
    do i=1,li
        code=cs(i)
        u=idx
        l=idx+lens(code)-1
        a=as(i)
        b=bs(i)
        v=vs(i)
        call gen_oh(code,a,b,v,xs(u:l),ys(u:l),zs(u:l),ws(u:l),lens(code))
        idx=idx+lens(code)
    end do
end subroutine ldfunc

LDFUNS
end module Lebdev
"""

temp="""
subroutine FUNC(xs,ys,zs,ws)
    implicit none
    real(8),intent(out)::xs(LDN),ys(LDN),zs(LDN),ws(LDN)
    real(8)::as(LEN),bs(LEN),vs(LEN)
    integer::cs(LEN)
    integer::li=LEN
    as=[&
    As
    bs=[&
    Bs
    vs=[&
    Vs
    cs=[&
    Cs
    call ldfunc(as,bs,vs,cs,li,xs,ys,zs,ws,LDN)
end subroutine FUNC
"""
tab='    '

def float_str(xs):
    res=''
    for i,x in enumerate(xs):
        res+=f'{x:<23.16E}'
        if i!=len(xs)-1:res+=','
        if (i+1)%5==0:res+=f'&\n{tab}'
    res=f'{res}]'
    return res.replace('E','D')

def intge_str(xs):
    res=''
    for i,x in enumerate(xs):
        res+=f'{x:<3.0f}'
        if i!=len(xs)-1:res+=','
        if (i+1)%20==0:res+=f'&\n{tab}'
    res=f'{res}]'
    return res


texts=[]

for i,func in enumerate(dir(lebedev)):
    if re.match(r'LD\d{4}',func) is None:continue
    paras:np.ndarray=eval(f'lebedev.{func}()')
    As=float_str(paras[:,0])
    Bs=float_str(paras[:,1])
    Vs=float_str(paras[:,2])
    Cs=intge_str(paras[:,3])
    ldn=int(func[2:])
    fstr=temp.replace('As',As)
    fstr=fstr.replace('Bs',Bs)
    fstr=fstr.replace('Vs',Vs)
    fstr=fstr.replace('Cs',Cs)
    fstr=fstr.replace('FUNC',func)
    fstr=fstr.replace('LEN',f'{paras.shape[0]}')
    fstr=fstr.replace('LDN',f'{ldn}')
    texts.append(fstr)

Path('lebdev.f90').write_text(text0.replace('LDFUNS',''.join(texts)),encoding='utf-8')