! 使用Fortran写一些比较耗时间的函数

module flib
    use iso_c_binding
contains

subroutine info() bind(c, name='info_')
    write(*,*)'Hello from Fortran'
end subroutine info

! 生成网格点坐标
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

! 计算波函数值
! 都使用fortran了，就不要用向量化计算了(numpy得用(╯▔皿▔)╯)，享受逐元素计算的快乐吧
subroutine gtf(alp,ngrid,grids,coord,l,m,n, wfn0, wfn1, level) bind(C, name="gtf_") ! 计算某些点处的高斯函数值，基函数
    ! 高斯指数，坐标数量，坐标值，坐标平方和，角动量分量，返回值
    ! 直接传入平方和防止重复计算，减少计算量
    use iso_c_binding
    implicit none
    integer(c_int),intent(in),value::ngrid
    real(c_double), intent(in) :: grids(3,ngrid)
    real(c_double),intent(in) :: coord(3) ! 原子坐标映射到基函数
    integer(c_int), intent(in), value :: l,m,n
    real(c_double), intent(in), value :: alp
    real(c_double), intent(inout) :: wfn0(ngrid)  ! 原始波函数值
    real(c_double), intent(inout) :: wfn1(3,ngrid)! 波函数一阶导
    integer(c_int), intent(in), value :: level ! 0表示不计算导数

    real(c_double)::pi = 3.1415926 ! 物理常量
    real(c_double)::Nm ! 归一化系数
    real(c_double)::r2 ! 坐标平方
    real(c_double)::fac ! 双阶乘
    real(c_double)::facs(0:2)
    integer(c_int)::ang ! 角动量
    real(c_double)::x,y,z
    real(c_double)::exv
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
        exv = exp(-alp*r2)
        wfn0(i) = Nm* x**l * y**m * z**n * exv
        if (level==0) continue ! 不计算导数
        !计算波函数梯度
        wfn1(1,i)=-2*alp*x*wfn0(i)
        wfn1(2,i)=-2*alp*y*wfn0(i)
        wfn1(3,i)=-2*alp*z*wfn0(i)

        if (l==1) wfn1(1,i)=wfn1(1,i) + y**m * z**n *exv*Nm
        if (m==1) wfn1(2,i)=wfn1(2,i) + x**l * z**n *exv*Nm
        if (n==1) wfn1(3,i)=wfn1(3,i) + x**l * y**m *exv*Nm

        if (l==2) wfn1(1,i)=wfn1(1,i) + x * y**m * z**n *exv*Nm*2
        if (m==2) wfn1(2,i)=wfn1(2,i) + x**l * y * z**n *exv*Nm*2
        if (n==2) wfn1(3,i)=wfn1(3,i) + x**l * y**m * z *exv*Nm*2
    end do
    
end subroutine gtf

! 计算某点处一个原子轨道波函数
subroutine cgf(cmax,nc,alps, coes,ngrid,grids,coord,l,m,n, wfn0, wfn1,level) bind(C, name="cgf_") ! 计算收缩波函数，原子轨道
    ! 收缩数量，高斯指数，收缩系数，角动量分量，点数量，点平方和，点坐标，波函数值
    use iso_c_binding
    implicit none
    integer(c_int), intent(in), value::nc,cmax ! 收缩数量
    real(c_double), intent(in) :: alps(cmax), coes(cmax)
    integer(c_int),intent(in),value::ngrid
    real(c_double), intent(in) :: grids(3,ngrid)
    real(c_double),intent(in) :: coord(3) ! 原子坐标映射到基函数
    integer(c_int), intent(in), value ::  l,m,n
    real(c_double), intent(inout) ::  wfn0(ngrid)
    real(c_double), intent(inout) ::  wfn1(3,ngrid)
    integer(c_int), intent(in), value :: level ! 0表示不计算导数
    real(c_double)::alp, coe, val0(ngrid), val1(3,ngrid)
    integer::i
    wfn0 = 0.0
    do i = 1, nc
        alp = alps(i)
        coe = coes(i)
        val0 = 0.0
        call gtf(alp, ngrid,grids,coord,l,m,n, val0, val1,level)
        wfn0 = wfn0 + coe*val0
        wfn1 = wfn1 + coe*val1
    end do
end subroutine cgf

! 计算所有基函数的原子轨道(向量)，分子轨道是原子轨道的线性组合，计算波函数及梯度
subroutine atoWfns(ngrid,grids,nmat,cords,cmax,ncgs,alpl,coel,lmns,wfn0,wfn1,level) bind(C, name="atoWfns_")
    implicit none
    integer(c_int), intent(in),value :: ngrid,nmat,cmax
    real(c_double), intent(in) :: grids(3,ngrid)
    real(c_double), intent(in) :: cords(3,nmat)
    real(c_double), intent(in) :: alpl(cmax,nmat)
    real(c_double), intent(in) :: coel(cmax,nmat)
    integer(c_int), intent(in) :: lmns(3,nmat) !角动量分量
    integer(c_int), intent(in) :: ncgs(nmat) ! 每一个原子轨道的收缩数量
    real(c_double), intent(inout)::wfn0(ngrid,nmat)
    real(c_double), intent(inout)::wfn1(3,ngrid,nmat)
    integer(c_int), intent(in), value :: level ! 0表示不计算导数

    integer(c_int)::l,m,n
    integer(c_int)::nc
    ! real(c_double) ::val0(ngrid), val1(3,ngrid)
    real(c_double) ::coord(3)
    real(c_double) ::alps(cmax),coes(cmax)
    integer::i

    ! write(*,*)'ngrid,nmat,cmax',ngrid,nmat,cmax

    do i=1,nmat
        l=lmns(1,i)
        m=lmns(2,i)
        n=lmns(3,i)
        coord=cords(:,i)
        alps=alpl(:,i)
        coes=coel(:,i)
        nc=ncgs(i)
        call cgf(cmax,nc,alps,coes,ngrid,grids,coord,l,m,n,wfn0(:,i),wfn1(:,:,i),level)
    end do
end subroutine atoWfns

! 计算分子电子密度及梯度，轨道的电子密度直接就是波函数的平方
subroutine molDens(ngrid,nmat,nobt,matC,wfns0,wfns1,dens0,dens1,level) bind(C,name="moldens_")
    use iso_c_binding
    integer(c_int),intent(in),value::ngrid
    integer(c_int),intent(in),value :: nmat ! 原子轨道数量，每一个原子轨道对应一个cgf
    integer(c_int),intent(in),value :: nobt ! 占据轨道数量，稀疏矩阵的列数
    real(c_double),intent(in) :: matC(nobt,nmat) ! 轨道系数矩阵
    real(c_double),intent(inout)::wfns0(ngrid,nmat) !波函数
    real(c_double),intent(inout)::wfns1(3,ngrid,nmat) !波函数梯度
    real(c_double),intent(inout)::dens0(ngrid)
    real(c_double),intent(inout)::dens1(3,ngrid)
    integer(c_int),intent(in),value :: level ! 0表示不计算导数
    real(c_double) ::wfn0(ngrid),wfn1(3,ngrid)
    integer :: obt,ato
    ! 提前算出所有原子轨道的波函数并存储起来，分子轨道的波函数只是原子轨道波函数的线性组合
    dens0=0.0
    dens1=0.0
    do obt=1,nobt !循环每一个分子轨道
        wfn0=0.0
        wfn1=0.0
        do ato=1,nmat ! 循环每一个原子轨道
            wfn0 = wfn0 + wfns0(:,ato)*matC(obt,ato)
            wfn1 = wfn1 + wfns1(:,:,ato)*matC(obt,ato)
        end do
        dens0=dens0+wfn0*wfn0
        dens1=dens1+wfn0*wfn1
    end do
end subroutine molDens

! 计算原子电子密度
subroutine atmDens(ngrid,nmat,matP,u,l,wfns,dens) bind(C,name="atmdens_")
    integer(c_int),intent(in),value::ngrid ! 网格数量
    integer(c_int),intent(in),value::nmat ! 矩阵大小
    integer(c_int),intent(in),value::u,l ! 该原子在矩阵中的上下界
    real(c_double),intent(in):: matP(nmat,nmat)
    real(c_double),intent(in)::wfns(ngrid,nmat) ! 将每个原子轨道的波函数存储下来
    real(c_double),intent(inout)::dens(ngrid)

    integer::i,j

    dens=0.0
    do i=u,l
        do j=1,nmat
            dens=dens+wfns(:,i)*wfns(:,j)*matP(j,i)
        end do
    end do
    
end subroutine atmDens

! 计算原子在分子格点的权重
subroutine a2mWeight(atm,nGrid,atmGrid,atmWeit,natm,atmPos,atmRad,atmDis,a2mGrid,a2mWeit,total) bind(C,name="a2mWeight_")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in),value :: atm ! 第多少个原子
    integer(c_int), intent(in),value :: nGrid ! 点的数量
    real(c_double), intent(in) :: atmGrid(3,nGrid)
    real(c_double), intent(in) :: atmWeit(nGrid)
    integer(c_int), intent(in),value :: natm ! 原子数量
    real(c_double), intent(in) :: atmPos(3,natm) ! 原子坐标
    real(c_double), intent(in) :: atmRad(natm) ! 原子半径
    real(c_double), intent(in) :: atmDis(natm,natm) ! 原子距离矩阵
    real(c_double), intent(out):: a2mGrid(3,nGrid) ! 修改后的坐标
    real(c_double), intent(out):: a2mWeit(nGrid) ! 原子权重
    integer(c_int),intent(out)::total

    real(c_double) :: S_u(natm,natm)
    real(c_double) :: pi(3),pj(3) ! 2个原子的坐标
    real(c_double) :: gp(3) ! 格点的坐标
    real(c_double) :: ri,rj ! 格点到原子的距离
    real(c_double) :: miu_ij,chi,nu_ij,a_ij,u_ij,rat
    real(c_double) :: wt(natm),weit
    integer::g,i,j
    ! do i=1,natm ! 打印原子坐标
    !     write(*,*)'atmPos',atmPos(:,i)
    ! end do
    total=0
    do g=1,nGrid
        gp=atmGrid(:,g)
        S_u=1.0
        do i=1,natm
            pi=atmPos(:,i)
            ri=sqrt(dot_product(pi-gp,pi-gp)) !两点之间的距离
            do j=1,natm
                if (i==j) cycle
                pj=atmPos(:,j)
                rj=sqrt(dot_product(pj-gp,pj-gp)) !两点之间的距离
                ! write(*,*)'gp',gp
                ! write(*,*)'ip',pi
                ! write(*,*)'jp',pj
                ! write(*,*)'ri,rj',ri,rj
                miu_ij=(ri-rj)/atmDis(j,i)
                ! write(*,*)'miu_ij',g,i,j,miu_ij
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
                ! write(*,*)'i,j,nu_ij',i,j,nu_ij
            end do
        end do
        ! write(*,*)'S_u',S_u
        wt=1.0
        do i=1,natm
            wt=wt*S_u(i,:)
        end do
        ! write(*,*)'wt',wt
        rat=wt(atm)/sum(wt)
        ! write(*,*)'atm,g,rat',atm,g,rat
        weit=atmWeit(g)*rat
        ! if (abs(weit)>1e-7) then
        if (.true.) then
            a2mWeit(g)=weit
            a2mGrid(:,g)=atmGrid(:,g)
            total=total+1
        end if
        
    end do
end subroutine a2mWeight

! 计算电子分布矩阵
subroutine eleMat(nmat,nobt,CM,SM,NM) bind(c, name="eleMat_")
    use iso_c_binding
    implicit none
    integer(c_int), intent(in),value :: nmat,nobt
    real(c_double), intent(in) ::  CM(nobt,nmat) ! 系数矩阵
    real(c_double), intent(in) ::  SM(nmat,nmat) ! 重叠矩阵
    real(c_double), intent(inout) ::  NM(nobt,nmat) ! 电子数量矩阵
    integer::i,j,b
    real(c_double)::val
    ! write(*,*)nmat,nobt
    ! write(*,*)'CM',CM(:,1)
    ! write(*,*)'SM',SM(:,1)
    NM=0.0
    do b=1,nmat
        do j=1,nobt
            val=0.0
            do i=1,nmat
                val=val+CM(j,b)*CM(j,i)*SM(b,i)
            end do
            NM(j,b)=val
        end do
    end do
        
end subroutine eleMat

subroutine nucPotential(ncord, cords, natm, nucs, xyzs, vals) bind(c, name="nucPotential_")!计算原子势
    integer(c_int), intent(in),value :: ncord !要计算的点的数量
    real(c_double), intent(out) ::  cords(3,ncord) !要计算的点的坐标
    integer(c_int), intent(in),value :: natm !原子数量
    integer(c_int), intent(in) :: nucs(natm) !原子核电荷
    real(c_double), intent(in) :: xyzs(3,natm) !原子坐标
    real(c_double), intent(out) :: vals(ncord) !要计算的点的势能值
    integer::i,j
    real(c_double)::dist !原子与格点之间的距离
    do i=1,ncord
        do j=1,natm
            dist=sqrt(dot_product(xyzs(:,j)-cords(:,i),xyzs(:,j)-cords(:,i))) !原子与格点之间的距离
            if (dist<1e-6) continue
            vals(i)=vals(i)+nucs(j)/dist
        end do
    end do
end subroutine nucPotential

subroutine elePotential(ncord,cords,ngrid, grids, weits, dens, vals) bind(c, name="elePotential_")!计算电子势
    integer(c_int), intent(in),value :: ncord !要计算的点的数量
    real(c_double), intent(in) ::  cords(3,ncord) !要计算的点的坐标
    integer(c_int), intent(in),value :: ngrid !dft格点的数量
    real(c_double), intent(in) ::  grids(3,ngrid) !dft格点
    real(c_double), intent(in) :: weits(ngrid) !格点权重
    real(c_double), intent(in) :: dens(ngrid) !电子密度
    
    real(c_double), intent(out) :: vals(ncord) !要计算的点的势能值
    integer::i,j
    real(c_double)::dists(ngrid) !坐标与格点之间的距离矩阵
    real(c_double)::wdens(ngrid)

    wdens=dens*weits !计算格点对应的电子密度x权重，减少计算量
    do i=1,ncord
        do j=1,ngrid
            dists(j)=sum((grids(:,j)-cords(:,i))**2)**0.5 ! 计算当前坐标点与所有格点之间的距离
        end do
        where (dists < 1e-6)
            dists = 1e-6
        end where
        vals(i)=sum(wdens/dists)
    end do
end subroutine elePotential

end module flib