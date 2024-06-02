! 使用Fortran写一些比较耗时间的函数

module flib
    use iso_c_binding
contains

subroutine add(x, y, res) bind(c, name='add_')
    use iso_c_binding
    implicit none
    real(c_float), intent(in), value :: x, y
    real(c_float), intent(out) ::  res
    res = x + y
end subroutine add

function sum2(a) result(b) bind(c, name='sum2_')
    use iso_c_binding
    implicit none

    real(c_float), intent(in)  :: a
    real(c_float)              :: b

    b = a + 2.0

end function sum2

subroutine double_array(x, N) bind(C, name="double_array_")
    use iso_c_binding
    implicit none

    integer(c_long), intent(in), value   :: N
    real(c_float), intent(inout)       :: x(N, N)

    x = exp(x)

end subroutine double_array

subroutine grid_pos(Nx, Ny, Nz, pos) bind(C, name="grid_pos_")
    use iso_c_binding
    implicit none
    integer(c_long), intent(in), value :: Nx, Ny, Nz
    real(c_float), intent(out) ::  pos(Nx*Ny*Nz, 3)
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
    integer(c_long), intent(in), value :: row, col
    real(c_float), intent(out) ::  pos(row, col)
    integer::i, j
    do i = 1, row
        do j = 1, col
            pos(i, j) = pos(i, j) + pos(i, j)
        end do
    end do
end subroutine same_array


! 都使用fortran了，就不要用向量化计算了(numpy得用(╯▔皿▔)╯)，享受逐元素计算的快乐吧
subroutine gtf(alp,ngrid,grids,coord,l,m,n, vals) bind(C, name="gtf_") ! 计算某些点处的高斯函数值，基函数
    ! 高斯指数，坐标数量，坐标值，坐标平方和，角动量分量，返回值
    ! 直接传入平方和防止重复计算，减少计算量
    use iso_c_binding
    implicit none
    integer(c_long),intent(in),value::ngrid
    real(c_float), intent(in) :: grids(3,ngrid)
    real(c_float),intent(in) :: coord(3) ! 原子坐标映射到基函数
    integer(c_long), intent(in), value :: l,m,n
    real(c_float), intent(in), value :: alp
    real(c_float), intent(inout) :: vals(ngrid)

    real(c_float)::pi = 3.1415926 ! 物理常量
    real(c_float)::Nm ! 归一化系数
    real(c_float)::r2 ! 坐标平方
    real(c_float)::fac ! 双阶乘
    real(c_float)::facs(0:2)
    integer(c_long)::ang ! 角动量
    real(c_float)::x,y,z
    real(c_float)::val
    integer::i

    facs=[1.,1.,3.]
    fac = facs(l)*facs(m)*facs(n) ! 数组越界不报错？
    ang = l + m + n
    Nm = (2.*alp/pi)**0.75*sqrt((4.*alp)**ang/fac) ! 计算归一化系数
    do i=1,ngrid
        
        x = grids(1,i)-coord(1)
        y = grids(2,i)-coord(2)
        z = grids(3,i)-coord(3)
        r2 = x**2 + y**2 + z**2
        val=x**l * y**m * z**n * exp(-alp*r2)*Nm
        vals(i) = val
        ! if ( i<=1 ) then
        !     write(*,*)'gtf',i,x,y,z,r2,Nm,val
        ! end if
        
    end do
    ! write(*,*)'gtf',sum(vals**2),fac,Nm,ngrid,alp
    
end subroutine gtf

! 计算某点处一个原子轨道波函数
subroutine cgf(cmax,nc,alps, coes,ngrid,grids,coord,l,m,n, wfn) bind(C, name="cgf_") ! 计算收缩波函数，原子轨道
    ! 收缩数量，高斯指数，收缩系数，角动量分量，点数量，点平方和，点坐标，波函数值
    use iso_c_binding
    implicit none
    integer(c_long), intent(in), value::nc,cmax ! 收缩数量
    real(c_float), intent(in) :: alps(cmax), coes(cmax)
    integer(c_long),intent(in),value::ngrid
    real(c_float), intent(in) :: grids(3,ngrid)
    real(c_float),intent(in) :: coord(3) ! 原子坐标映射到基函数
    integer(c_long), intent(in), value ::  l,m,n
    real(c_float), intent(inout) ::  wfn(ngrid)
    real(c_float)::alp, coe, vals(ngrid)
    integer::i

    wfn = 0.0
    do i = 1, nc
        alp = alps(i)
        coe = coes(i)
        vals = 0.0
        call gtf(alp, ngrid,grids,coord,l,m,n, vals)
        wfn = wfn + coe*vals
    end do
    ! write(*,*)'cgf',sum(wfn**2)
    ! write(*,*)'cgf done'
end subroutine cgf



! 计算某点处所有原子轨道波函数加和
subroutine obtWfn(ngrid,grids,nmat,cords,ncgs,cmax,oalps,ocoes,coefs,lmns,wfn) bind(C,name="obtwfn_")
    use iso_c_binding
    integer(c_long),intent(in),value::ngrid
    real(c_float), intent(in) :: grids(3,ngrid)
    integer(c_long),intent(in),value :: nmat ! 原子轨道数量，每一个原子轨道对应一个cgf
    integer(c_long),intent(in) :: ncgs(nmat) ! 每一个原子轨道的收缩数量
    real(c_float),intent(in) :: cords(3,nmat) ! 原子坐标映射到基函数
    integer(c_long),intent(in),value :: cmax !最大收缩数量
    real(c_float),intent(in) :: oalps(cmax,nmat)
    real(c_float),intent(in) :: ocoes(cmax,nmat)
    real(c_float),intent(in) :: coefs(nmat)
    integer(c_long),intent(in) :: lmns(3,nmat) !角动量分量
    real(c_float),intent(inout)::wfn(ngrid)
    
    real(c_float)::alps(cmax) ! 高斯指数
    real(c_float)::coes(cmax) ! 收缩系数
    integer(c_long)::nc ! 每个原子轨道收缩的数量
    integer(c_long)::l,m,n
    real(c_float)::val(ngrid)
    real(c_float)::coord(3)

    integer::i

    wfn=0.0
    do i=1,nmat
        nc=ncgs(i)
        ! write(*,*)i,nc
        alps=oalps(:,i)
        coes=ocoes(:,i)
        ! write(*,*)'alps',alps
        l=lmns(1,i)
        m=lmns(2,i)
        n=lmns(3,i)
        coord=cords(:,i)
        val=0.0
        call cgf(cmax,nc,alps,coes,ngrid,grids,coord,l,m,n,val)
        ! write(*,*)'obtWfn_i',x,y,z,val
        wfn = wfn + coefs(i)*val
    end do
    ! write(*,*)'obtWfn',sum(wfn**2),coefs
end subroutine obtWfn

! 计算某一点处的分子电子密度，轨道的电子密度直接就是波函数的平方
subroutine molDens(ngrid,grids,nmat,cords,nobt,CM,ncgs,cmax,oalps,ocoes,lmns,dens) bind(C,name="moldens_")
    use iso_c_binding
    integer(c_long),intent(in),value::ngrid
    real(c_float), intent(in) :: grids(3,ngrid)
    integer(c_long),intent(in),value :: nmat ! 原子轨道数量，每一个原子轨道对应一个cgf
    real(c_float),intent(in) :: cords(3,nmat) ! 原子坐标映射到基函数
    integer(c_long),intent(in),value :: nobt ! 占据轨道数量，稀疏矩阵的列数
    real(c_float),intent(in) :: CM(nobt,nmat) ! 轨道系数矩阵
    integer(c_long),intent(in) :: ncgs(nmat) ! 每一个原子轨道的收缩数量
    integer(c_long),intent(in),value :: cmax !最大收缩数量
    real(c_float),intent(in) :: oalps(cmax,nmat)
    real(c_float),intent(in) :: ocoes(cmax,nmat)
    integer(c_long),intent(in) :: lmns(3,nmat) !角动量分量
    real(c_float),intent(inout)::dens(ngrid)

    real(c_float)::coefs(nmat)
    real(c_float) ::wfn(ngrid)
    integer :: i
    ! write(*,*)'molDen',nmat,nobt,cmax
    ! do i=1,nmat
    !     write(*,*)'CMi',i,CM(:,i)
    ! end do
    dens=0.0
    do i=1,nobt
        coefs = CM(i,:)
        wfn = 0.0
        call obtwfn(ngrid,grids,nmat,cords,ncgs,cmax,oalps,ocoes,coefs,lmns,wfn)
        dens = dens + wfn**2
    end do
end subroutine molDens


subroutine a2mWeight(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis,a2mGrid,a2mWeit) bind(C,name="a2mWeight_")
    use iso_c_binding
    implicit none
    integer(c_long), intent(in),value :: atm ! 第多少个原子
    integer(c_long), intent(in),value :: nGrid ! 点的数量
    real(c_float), intent(in) :: atmGrid(3,nGrid)
    real(c_float), intent(in) :: atmWeit(nGrid)
    integer(c_long), intent(in),value :: natm ! 原子数量
    real(c_float), intent(in) :: atmPos(3,natm) ! 原子坐标
    real(c_float), intent(in) :: atmRad(natm) ! 原子半径
    real(c_float), intent(in) :: atmDis(natm,natm) ! 原子距离矩阵
    real(c_float), intent(out):: a2mGrid(3,nGrid) ! 修改后的坐标
    real(c_float), intent(out):: a2mWeit(nGrid) ! 原子权重

    real(c_float) :: S_u(natm,natm)
    real(c_float) :: pi(3),pj(3) ! 量个原子的坐标
    real(c_float) :: gp(3) ! 格点的坐标
    real(c_float) :: ri,rj ! 格点到原子的距离
    real(c_float) :: miu_ij,chi,nu_ij,a_ij,u_ij,rat
    real(c_float) :: wt(natm)
    integer::g,i,j
    ! do i=1,natm
    !     write(*,*)atmPos(:,i)
    ! end do
    do g=1,nGrid
        gp=atmGrid(:,g)
        S_u=1.0
        do i=1,natm
            pi=atmPos(:,i)

            ri=sqrt(dot_product(pi-gp,pi-gp)) !量向量之间的距离
            do j=1,natm
                if (i==j) cycle
                pj=atmPos(:,j)
                rj=sqrt(dot_product(pj-gp,pj-gp)) !量向量之间的距离

                miu_ij=(ri-rj)/atmDis(j,i)
                chi=atmRad(i)/atmRad(j) !两原子半径之比
                if (abs(chi-1)<1e-6) then
                    nu_ij=miu_ij
                else
                    u_ij=(chi-1)/(chi+1)
                    a_ij=u_ij/(u_ij**2-1)
                    if (a_ij>0.5) a_ij=0.5
                    if (a_ij<-0.5) a_ij=-0.5
                    nu_ij=miu_ij+a_ij*(1-miu_ij**2)
                end if
                
                nu_ij=1.5*nu_ij-0.5*nu_ij**3
                nu_ij=1.5*nu_ij-0.5*nu_ij**3
                nu_ij=1.5*nu_ij-0.5*nu_ij**3

                S_u(j,i)=0.5*(1-nu_ij)
                
            end do
        end do
        wt=1.0
        do i=1,natm
            wt=wt*S_u(i,:)
        end do

        rat=wt(atm)/sum(wt)
        a2mWeit(g)=atmWeit(g)*rat
        a2mGrid(:,g)=atmGrid(:,g)
    end do
end subroutine a2mWeight

end module flib