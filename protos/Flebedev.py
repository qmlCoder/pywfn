"""
将python的lebedev转为Fortran的
"""
import sys
sys.path.append('d:/code/pywfn')
from pywfn.data import lebedev
import re
from pathlib import Path
import numpy as np

text0="""
module Lebdev
    implicit none
    integer,parameter::cnum(6)=[6,12,8,24,24,48] !每个code对应的点数
contains
subroutine gen_oh(code,a,b,v,grids,weits,n)
    integer, intent(in) :: code
    integer, intent(in) :: n !生成格点数量
    real(8), intent(inout) :: a,b
    real(8), intent(in) :: v
    real(8), intent(out):: grids(3,n),weits(n)
    real(8)::z=0.0d0,c
    select case(code)
    case(1)
        a=1.0d0
        grids(1,:)=[ a,-a, z, z, z, z]
        grids(2,:)=[ z, z, a,-a, z, z]
        grids(3,:)=[ z, z, z, z, a,-a]
        weits     =[ v, v, v, v, v, v]
    case(2)
        a=sqrt(1.0/2.0)
        grids(1,:)=[ z, z, z, z, a,-a, a,-a, a,-a, a,-a]
        grids(2,:)=[ a,-a, a,-a, z, z, z, z, a, a,-a,-a]
        grids(3,:)=[ a, a,-a,-a, a, a,-a,-a, z, z, z, z]
        weits     =[ v, v, v, v, v, v, v, v, v, v, v, v]
    case(3)
        a=sqrt(1.0/3.0)
        grids(1,:)=[ a,-a, a,-a, a,-a, a,-a]
        grids(2,:)=[ a, a,-a,-a, a, a,-a,-a]
        grids(3,:)=[ a, a, a, a,-a,-a,-a,-a]
        weits     =[ v, v, v, v, v, v, v, v]
    case(4)
        b = sqrt(1.0-2*a**2)
        grids(1,:)=[ a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, b,-b, b,-b, b,-b, b,-b]
        grids(2,:)=[ a, a,-a,-a, a, a,-a,-a, b, b,-b,-b, b, b,-b,-b, a, a,-a,-a, a, a,-a,-a]
        grids(3,:)=[ b, b, b, b,-b,-b,-b,-b, a, a, a, a,-a,-a,-a,-a, a, a, a, a,-a,-a,-a,-a]
        weits     =[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    case(5)
        b=sqrt(1.0-a**2)
        grids(1,:)=[ a,-a, a,-a, b,-b, b,-b, a,-a, a,-a, b,-b, b,-b, z, z, z, z, z, z, z, z]
        grids(2,:)=[ b, b,-b,-b, a, a,-a,-a, z, z, z, z, z, z, z, z, a,-a, a,-a, b,-b, b,-b]
        grids(3,:)=[ z, z, z, z, z, z, z, z, b, b,-b,-b, a, a,-a,-a, b, b,-b,-b, a, a,-a,-a]
        weits     =[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    case (6)
        c=sqrt(1.0-a**2-b**2)
        grids(1,:)=[ a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, a,-a, b,-b, b,-b, b,-b, b,-b,&
                     b,-b, b,-b, b,-b, b,-b, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c, c,-c]
        grids(2,:)=[ b, b,-b,-b, b, b,-b,-b, c, c,-c,-c, c, c,-c,-c, a, a,-a,-a, a, a,-a,-a,&
                     c, c,-c,-c, c, c,-c,-c, a, a,-a,-a, a, a,-a,-a, b, b,-b,-b, b, b,-b,-b]
        grids(3,:)=[ c, c, c, c,-c,-c,-c,-c, b, b, b, b,-b,-b,-b,-b, c, c, c, c,-c,-c,-c,-c,&
                     a, a, a, a,-a,-a,-a,-a, b, b, b, b,-b,-b,-b,-b, a, a, a, a,-a,-a,-a,-a]
        weits     =[ v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v,&
                     v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v, v]
    end select
end subroutine gen_oh

subroutine ldfunc(as,bs,vs,cs,ni,grids,weits,no)
    integer,intent(in)::ni,no
    real(8), intent(in) :: as(ni),bs(ni),vs(ni)
    integer, intent(in) :: cs(ni) !code
    real(8), intent(out) :: grids(3,no),weits(no)
    real(8)::a,b,v
    integer::i
    integer::up,lo !索引的上下界
    integer::code
    
    integer::idx=1
    do i=1,ni
        code=cs(i)
        up=idx
        lo=idx+cnum(code)-1
        a=as(i)
        b=bs(i)
        v=vs(i)
        call gen_oh(code,a,b,v,grids(:,up:lo),weits(up:lo),cnum(code))
        idx=idx+cnum(code)
    end do
    if (idx-1/=no) then
        write(*,*)"lebdev 数量不对" idx-1,no
        stop
    end if
end subroutine ldfunc

LDFUNS
end module Lebdev
"""

temp="""
subroutine FUNC(grids,weits)
    implicit none
    real(8),intent(out)::grids(3,LDN),weits(LDN)
    real(8)::as(LEN),bs(LEN),vs(LEN)
    integer::cs(LEN)
    integer::ni=LEN
    as=[&
    As
    bs=[&
    Bs
    vs=[&
    Vs
    cs=[&
    Cs
    call ldfunc(as,bs,vs,cs,ni,grids,weits,LDN)
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