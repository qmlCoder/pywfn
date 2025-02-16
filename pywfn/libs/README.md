```fortran
    !   ! dxx
    !   select case (l)
    !   case (0)
    !     wfn2(1, 1, i) = wfn0(i)*2*alp*(2*x**2*alp - 1)
    !   case (1)
    !     wfn2(1, 1, i) = wfn0(i)*(-4*alp + 2*alp*(2*x**2*alp - 1))
    !   case (2)
    !     wfn2(1, 1, i) = wfn0(i)*(-8*alp + 2*alp*(2*x**2*alp - 1)) + 2*y**m*z**n*exv*Nm
    !   end select
    !   ! dyy
    !   select case (m)
    !   case (0)
    !     wfn2(2, 2, i) = wfn0(i)*2*alp*(2*y**2*alp - 1)
    !   case (1)
    !     wfn2(2, 2, i) = wfn0(i)*(-4*alp + 2*alp*(2*y**2*alp - 1))
    !   case (2)
    !     wfn2(2, 2, i) = wfn0(i)*(-8*alp + 2*alp*(2*y**2*alp - 1)) + 2*x**l*z**n*exv*Nm
    !   end select
    !   ! dzz
    !   select case (n)
    !   case (0)
    !     wfn2(3, 3, i) = wfn0(i)*2*alp*(2*z**2*alp - 1)
    !   case (1)
    !     wfn2(3, 3, i) = wfn0(i)*(-4*alp + 2*alp*(2*z**2*alp - 1))
    !   case (2)
    !     wfn2(3, 3, i) = wfn0(i)*(-8*alp + 2*alp*(2*z**2*alp - 1)) + 2*x**l*y**m*exv*Nm
    !   end select
    !   ! dxy
    !   if (l == 0 .and. m == 0) then
    !     wfn2(2, 1, i) = z**n*4*alp**2*x*y*exv*Nm
    !   else if (l == 0) then
    !     wfn2(2, 1, i) = y**(m - 1)*z**n*(4*alp**2*x*y**2 - 2*alp*m*x)*exv*Nm
    !   else if (m == 0) then
    !     wfn2(2, 1, i) = x**(l - 1)*z**n*(4*alp**2*x**2*y - 2*alp*l*y)*exv*Nm
    !   else
    !     wfn2(2, 1, i) = x**(l - 1)*y**(m - 1)*z**n*(4*alp**2*x**2*y**2 - 2*alp*l*y**2 - 2*alp*m*x**2 + l*m)*exv*Nm
    !   end if
    !   ! dxz
    !   if (l == 0 .and. n == 0) then
    !     wfn2(3, 1, i) = y**m*4*alp**2*x*z*exv*Nm
    !   else if (l == 0) then
    !     wfn2(3, 1, i) = y**m*z**(n - 1)*(4*alp**2*x*z**2 - 2*alp*n*x)*exv*Nm
    !   else if (n == 0) then
    !     wfn2(3, 1, i) = x**(l - 1)*y**m*(4*alp**2*x**2*z - 2*alp*l*z)*exv*Nm
    !   else
    !     wfn2(3, 1, i) = x**(l - 1)*y**m*z**(n - 1)*(4*alp**2*x**2*z**2 - 2*alp*l*z**2 - 2*alp*n*x**2 + l*n)*exv*Nm
    !   end if
    !   ! dyz
    !   if (m == 0 .and. n == 0) then
    !     wfn2(3, 2, i) = x**l*4*alp**2*y*z*exv*Nm
    !   else if (m == 0) then
    !     wfn2(3, 2, i) = x**l*z**(n - 1)*(4*alp**2*y*z**2 - 2*alp*n*y)*exv*Nm
    !   else if (n == 0) then
    !     wfn2(3, 2, i) = x**l*y**(m - 1)*(4*alp**2*y**2*z - 2*alp*m*z)*exv*Nm
    !   else
    !     wfn2(3, 2, i) = x**l*y**(m - 1)*z**(n - 1)*(4*alp**2*y**2*z**2 - 2*alp*m*z**2 - 2*alp*n*y**2 + m*n)*exv*Nm
    !   end if
    !   wfn2(1, 2, i) = wfn2(2, 1, i)
    !   wfn2(1, 3, i) = wfn2(3, 1, i)
    !   wfn2(2, 3, i) = wfn2(3, 2, i)

```

## 编译命令
```shell
gfortran -shared -ffree-form -ffree-line-length-none -fopenmp data.f90 flib.f90 -o flib.dll
```

./gfile/ch4_6d.fch