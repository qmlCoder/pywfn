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


! 都使用fortran了，就不要用向量化计算了(numpy得用(╯▔皿▔)╯)，享受逐元素计算的快乐吧
subroutine gtf(alp, x,y,z,l,m,n, val) bind(C, name="gtf_") ! 计算某些点处的高斯函数值，基函数
    ! 高斯指数，坐标数量，坐标值，坐标平方和，角动量分量，返回值
    ! 直接传入平方和防止重复计算，减少计算量
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value :: x,y,z
    integer(c_int), intent(in), value :: l,m,n
    real(c_double), intent(in), value :: alp
    real(c_double), intent(out) :: val

    real(c_double)::pi = 3.1415926 ! 物理常量
    real(c_double)::N_ ! 归一化系数
    real(c_double)::r2 ! 归一化系数
    integer(c_int)::fl, fm, fn, fac ! 双阶乘
    integer(c_int)::ang ! 角动量


    call fac2(2*l - 1, fl)
    call fac2(2*m - 1, fm)
    call fac2(2*n - 1, fn)
    fac = fl*fm*fn
    ang = l + m + n
    N_ = (2*alp/pi)**(3/4)*sqrt((4*alp)**ang/fac) ! 计算归一化系数
    r2=x**2+y**2+z**2
    val = x**l*y**m*z**n*exp(-alp*r2)*N_
    
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

subroutine a2mWeight(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis,a2mGrid,a2mWeit) bind(C,name="a2mWeight_")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in) :: atm ! 第多少个原子
    integer(c_int), intent(in) :: nGrid ! 点的数量
    real(c_double), intent(in) :: atmGrid(nGrid,3)
    real(c_double), intent(in) :: atmWeit(nGrid)
    integer(c_int), intent(in) :: natm ! 原子数量
    real(c_double), intent(in) :: atmPos(natm,3) ! 原子坐标
    real(c_double), intent(in) :: atmRad(natm) ! 原子半径
    real(c_double), intent(in) :: atmDis(natm,natm) ! 原子距离矩阵
    real(c_double), intent(out):: a2mGrid(nGrid,3) ! 修改后的坐标
    real(c_double), intent(out):: a2mWeit(nGrid) ! 原子权重

    real(c_double) :: S_u(natm,natm)
    real(c_double) :: pi(3),pj(3) ! 量个原子的坐标
    real(c_double) :: pg(3) ! 格点的坐标
    real(c_double) :: ri,rj ! 格点到原子的距离
    real(c_double) :: miu_ij,chi,nu_ij,a_ij,u_ij,rat
    real(c_double) :: wt(natm)
    integer::g,i,j
    do g=1,nGrid
        pg=atmGrid(g,:)
        S_u=0.0
        do i=1,natm
            pi=atmPos(i,:)
            ri=sqrt(dot_product(pi-pg,pi-pg)) !量向量之间的距离
            do j=1,natm
                pj=atmPos(j,:)
                rj=sqrt(dot_product(pj-pg,pj-pg)) !量向量之间的距离
                miu_ij=(ri-rj)/atmDis(i,j)
                chi=atmRad(i)/atmRad(j) !两原子半径之比
                if (abs(chi-1)<1e-6) then
                    nu_ij=miu_ij
                else
                    u_ij=(chi-1)/(chi+1)
                    a_ij=nu_ij/(u_ij**2-1)
                    if (a_ij>0.5) a_ij=0.5
                    if (a_ij<-0.5) a_ij=-0.5
                    nu_ij=miu_ij+a_ij*(1-miu_ij**2)
                
                nu_ij=1.5*nu_ij-0.5*nu_ij**3
                nu_ij=1.5*nu_ij-0.5*nu_ij**3
                nu_ij=1.5*nu_ij-0.5*nu_ij**3

                S_u(i,j)=0.5*(1-nu_ij)
                end if
            end do
        end do
        wt=1.0
        do i=1,natm
            wt=wt*S_u(i,:)
        rat=wt(atm)/sum(wt)
        a2mWeit(g)=atmWeit(g)*rat ! 修改权重
        a2mGrid(g,:)=atmGrid(g,:) ! 赋值坐标
        end do
    end do

end subroutine a2mWeight