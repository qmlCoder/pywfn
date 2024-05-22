! 使用Fortran写一些比较耗时间的函数

subroutine add(x, y, res) bind(c, name='add_')
    use iso_c_binding
    implicit none
    real(c_double), intent(in), value :: x, y
    real(c_double), intent(out) ::  res
    res = x + y
end subroutine add

function sum2(a) result(b) bind(c, name='sum2_')
    use iso_c_binding
    implicit none

    real(c_double), intent(in)  :: a
    real(c_double)              :: b

    b = a + 2.d0

end function sum2

subroutine double_array(x, N) bind(C, name="double_array_")
    use iso_c_binding
    implicit none

    integer(c_int), intent(in), value   :: N
    real(c_double), intent(inout)       :: x(N, N)

    x = exp(x)

end subroutine double_array

subroutine grid_pos(Nx, Ny, Nz, pos) bind(C, name="grid_pos_")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: Nx, Ny, Nz
    real(c_double), intent(out) ::  pos(Nx*Ny*Nz, 3)
    integer::i, j, k, l
    l = 1
    do i = 0, Nx - 1
        do j = 0, Ny - 1
            do k = 0, Nz - 1
                pos(l, :) = [i, j, k]
                l = l + 1
            end do
        end do
    end do

end subroutine grid_pos

subroutine same_array(row, col, pos) bind(C, name="same_array_")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: row, col
    real(c_double), intent(out) ::  pos(row, col)
    integer::i, j
    do i = 1, row
        do j = 1, col
            pos(i, j) = pos(i, j) + pos(i, j)
        end do
    end do
end subroutine same_array

subroutine fac2(num, res) bind(C, name="fac2_") !计算双阶乘
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: num
    integer(c_int), intent(out) :: res
    integer::i
    res = 1
    if (num > 1) then
        do i = 1, num, 2
            res = res*i
        end do
    end if
end subroutine fac2

subroutine gtf(alp, np, pos, r2, lmn, val) bind(C, name="gtf_") ! 计算某些点处的高斯函数值，基函数
    ! 高斯指数，坐标数量，坐标值，坐标平方和，角动量分量，返回值
    ! 直接传入平方和防止重复计算，减少计算量
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: np
    integer(c_int), intent(in) :: lmn(3)
    real(c_double), intent(in), value :: alp
    real(c_double), intent(in) :: pos(np, 3), r2(np)
    real(c_double), intent(out) :: val(np)
    real(c_double)::xs(np), ys(np), zs(np)
    integer(c_int)::l, m, n ! 角动量分量
    real(c_double)::pi = 3.1415926 ! 物理常量
    real(c_double)::N_ ! 归一化系数
    integer(c_int)::fl, fm, fn, fac ! 双阶乘
    integer(c_int)::ang ! 角动量

    l = lmn(1)
    m = lmn(2)
    n = lmn(3)
    call fac2(2*l - 1, fl)
    call fac2(2*m - 1, fm)
    call fac2(2*n - 1, fn)
    fac = fl*fm*fn
    ang = l + m + n
    N_ = (2*alp/pi)**(3/4)*sqrt((4*alp)**ang/fac) ! 计算归一化系数
    xs = pos(:, 1)
    ys = pos(:, 2)
    zs = pos(:, 3)
    write(*,*)size(val)
    val = xs**l*ys**m*zs**n*exp(-alp*r2)*N_
    
end subroutine gtf

subroutine cgf(nc, alps, coes, lmn, np, r2, pos, wfn) bind(C, name="cgf_") ! 计算收缩波函数，原子轨道
    ! 收缩数量，高斯指数，收缩系数，角动量分量，点数量，点平方和，点坐标，波函数值
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value::nc ! 收缩数量
    real(c_double), intent(in) :: alps(nc), coes(nc)
    integer(c_int), intent(in) ::  lmn(nc)
    integer(c_int), intent(in), value :: np
    real(c_double), intent(in) :: pos(np, 3), r2(np)
    real(c_double), intent(out) ::  wfn(nc)
    real(c_double)::alp, coe, val(np)
    integer::i
    val=1.0
    do i = 1, nc
        alp = alps(i)
        coe = coes(i)
        call gtf(alp, np, pos, r2, lmn, val)
        wfn = wfn + coe*val
    end do
end subroutine cgf
