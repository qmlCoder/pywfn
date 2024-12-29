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
    real(c_double), intent(out) ::  pos(3, Nx*Ny*Nz)
    integer::i, j, k, l
    l = 1
    !$omp parallel default(none), private(i,j,k,l), shared(pos,Nx,Ny,Nz)
    !$omp do
    do i = 1, Nx
        do j = 1, Ny
            do k = 1, Nz
                l = (i-1)*Ny*Nz + (j-1)*Nz + k
                pos(:, l) = [real(i-1, c_double), real(j-1, c_double), real(k-1, c_double)]
            end do
        end do
    end do
    !$omp end do
    !$omp end parallel
end subroutine grid_pos

! 计算波函数值
subroutine gtf(alp,ngrid,grids,coord,l,m,n,level, wfn0, wfn1, wfn2) bind(C, name="gtf_") ! 计算某些点处的高斯函数值，基函数
    ! 高斯指数，坐标数量，坐标值，坐标平方和，角动量分量，返回值
    ! 直接传入平方和防止重复计算，减少计算量
    use iso_c_binding
    implicit none
    integer(c_int),intent(in),value::ngrid
    real(c_double), intent(in) :: grids(3,ngrid)
    real(c_double),intent(in) :: coord(3) ! 原子坐标映射到基函数
    integer(c_int), intent(in), value :: l,m,n
    real(c_double), intent(in), value :: alp
    real(c_double), intent(inout) :: wfn0(ngrid) ! 原始波函数值
    real(c_double), intent(inout) :: wfn1(3,ngrid)! 波函数一阶导
    real(c_double), intent(inout) :: wfn2(3,3,ngrid) ! 波函数二阶导
    integer(c_int), intent(in), value :: level ! 0：不计算倒数，1：计算一介倒数，2：计算二阶导数

    ! real(c_double)::pi = 3.1415926 ! 物理常量
    real(c_double)::pi = 4.*atan(1.) ! 物理常量
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
    wfn0=0.0
    wfn1=0.0
    wfn2=0.0
    !$omp parallel do default(none), private(i,x,y,z,r2,exv), shared(ngrid,grids,coord,l,m,n,alp,wfn0,wfn1,wfn2,level,Nm)
    do i=1,ngrid
        x = grids(1,i)-coord(1)
        y = grids(2,i)-coord(2)
        z = grids(3,i)-coord(3)
        r2 = x**2 + y**2 + z**2
        exv = exp(-alp*r2)
        wfn0(i) = Nm* x**l * y**m * z**n * exv
        if (level==0) cycle ! 不计算导数
        ! 计算波函数梯度
        wfn1(1,i)=-2*alp*x*wfn0(i)
        wfn1(2,i)=-2*alp*y*wfn0(i)
        wfn1(3,i)=-2*alp*z*wfn0(i)

        if (l==1) wfn1(1,i)=wfn1(1,i) + y**m * z**n *exv*Nm
        if (m==1) wfn1(2,i)=wfn1(2,i) + x**l * z**n *exv*Nm
        if (n==1) wfn1(3,i)=wfn1(3,i) + x**l * y**m *exv*Nm

        if (l==2) wfn1(1,i)=wfn1(1,i) + x * y**m * z**n *exv*Nm*2
        if (m==2) wfn1(2,i)=wfn1(2,i) + x**l * y * z**n *exv*Nm*2
        if (n==2) wfn1(3,i)=wfn1(3,i) + x**l * y**m * z *exv*Nm*2

        if (level==1) cycle ! 不计算二阶导
        ! 计算波函数二阶导数
        wfn2(1,1,i) = Nm*x**(l - 2)*y**m*z**n*(2*alp*x**2*(2*alp*x**2 - 2*l - 1) + l*(l - 1))               *exp(-alp*(x**2 + y**2 + z**2))
        wfn2(1,2,i) = Nm*x**(l - 1)*y**(m - 1)*z**n*(4*alp**2*x**2*y**2 - 2*alp*l*y**2 - 2*alp*m*x**2 + l*m)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(1,3,i) = Nm*x**(l - 1)*y**m*z**(n - 1)*(4*alp**2*x**2*z**2 - 2*alp*l*z**2 - 2*alp*n*x**2 + l*n)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(2,1,i) = Nm*x**(l - 1)*y**(m - 1)*z**n*(4*alp**2*x**2*y**2 - 2*alp*l*y**2 - 2*alp*m*x**2 + l*m)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(2,2,i) = Nm*x**l*y**(m - 2)*z**n*(2*alp*y**2*(2*alp*y**2 - 2*m - 1) + m*(m - 1))               *exp(-alp*(x**2 + y**2 + z**2))
        wfn2(2,3,i) = Nm*x**l*y**(m - 1)*z**(n - 1)*(4*alp**2*y**2*z**2 - 2*alp*m*z**2 - 2*alp*n*y**2 + m*n)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(3,1,i) = Nm*x**(l - 1)*y**m*z**(n - 1)*(4*alp**2*x**2*z**2 - 2*alp*l*z**2 - 2*alp*n*x**2 + l*n)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(3,2,i) = Nm*x**l*y**(m - 1)*z**(n - 1)*(4*alp**2*y**2*z**2 - 2*alp*m*z**2 - 2*alp*n*y**2 + m*n)*exp(-alp*(x**2 + y**2 + z**2))
        wfn2(3,3,i) = Nm*x**l*y**m*z**(n - 2)*(2*alp*z**2*(2*alp*z**2 - 2*n - 1) + n*(n - 1))               *exp(-alp*(x**2 + y**2 + z**2))
    end do
    !$omp end parallel do
    
end subroutine gtf

! 计算某点处一个原子轨道波函数，收缩基函数
subroutine cgf(cmax,nc,alps, coes,ngrid,grids,coord,l,m,n, level, wfn0, wfn1, wfn2) bind(C, name="cgf_") ! 计算收缩波函数，原子轨道
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
    real(c_double), intent(inout) ::  wfn2(3,3,ngrid)
    integer(c_int), intent(in), value :: level ! 0表示不计算导数
    real(c_double)::alp, coe, val0(ngrid), val1(3,ngrid),val2(3,3,ngrid)
    integer::i
    wfn0 = 0.0
    do i = 1, nc
        alp = alps(i)
        coe = coes(i)
        val0 = 0.0
        call gtf(alp, ngrid,grids,coord,l,m,n, level ,val0, val1, val2)
        wfn0 = wfn0 + coe*val0
        wfn1 = wfn1 + coe*val1
        wfn2 = wfn2 + coe*val2
    end do
end subroutine cgf

! 计算所有基函数的原子轨道(向量)，分子轨道是原子轨道的线性组合，计算波函数及梯度
subroutine atoWfns(ngrid,grids,nmat,cords,cmax,ncgs,alpl,coel,lmns,level,wfn0,wfn1,wfn2) bind(C, name="atoWfns_")
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
    real(c_double), intent(inout)::wfn2(3,3,ngrid,nmat)
    integer(c_int), intent(in), value :: level ! 0表示不计算导数

    integer(c_int)::l,m,n
    integer(c_int)::nc
    ! real(c_double) ::val0(ngrid), val1(3,ngrid)
    real(c_double) ::coord(3)
    real(c_double) ::alps(cmax),coes(cmax)
    integer::i

    ! write(*,*)'ngrid,nmat,cmax',ngrid,nmat,cmax
    wfn0=0.0
    wfn1=0.0
    wfn2=0.0
    do i=1,nmat
        l=lmns(1,i)
        m=lmns(2,i)
        n=lmns(3,i)
        coord=cords(:,i)
        alps=alpl(:,i)
        coes=coel(:,i)
        nc=ncgs(i)
        call cgf(cmax,nc,alps,coes,ngrid,grids,coord,l,m,n,level,wfn0(:,i),wfn1(:,:,i),wfn2(:,:,:,i))
    end do
end subroutine atoWfns

! 计算分子电子密度及梯度，轨道的电子密度直接就是波函数的平方
subroutine molDens(ngrid,nmat,nobt,matC,wfns0,wfns1,wfns2,level,dens0,dens1,dens2) bind(C,name="moldens_")
    use iso_c_binding
    integer(c_int),intent(in),value::ngrid
    integer(c_int),intent(in),value :: nmat ! 原子轨道数量，每一个原子轨道对应一个cgf
    integer(c_int),intent(in),value :: nobt ! 占据轨道数量，稀疏矩阵的列数
    real(c_double),intent(in) :: matC(nobt,nmat) ! 轨道系数矩阵
    real(c_double),intent(inout)::wfns0(ngrid,nmat) !波函数
    real(c_double),intent(inout)::wfns1(3,ngrid,nmat) !波函数梯度
    real(c_double),intent(inout)::wfns2(3,3,ngrid,nmat) !波函数二阶导

    real(c_double),intent(inout)::dens0(ngrid)
    real(c_double),intent(inout)::dens1(3,ngrid)
    real(c_double),intent(inout)::dens2(3,3,ngrid)
    integer(c_int),intent(in),value :: level ! 0表示不计算导数
    real(c_double) ::wfn0(ngrid,nobt),wfn1(3,ngrid,nobt),wfn2(3,3,ngrid,nobt)
    real(c_double) ::den0(ngrid,nobt),den1(3,ngrid,nobt),den2(3,3,ngrid,nobt)
    integer :: obt,ato,i,j
    ! 提前算出所有原子轨道的波函数并存储起来，分子轨道的波函数只是原子轨道波函数的线性组合
    dens0=0.0
    dens1=0.0
    dens2=0.0

    wfn0=0.0
    wfn1=0.0
    wfn2=0.0
    den0=0.0
    den1=0.0
    den2=0.0
    do obt=1,nobt !循环每一个分子轨道
        !$omp parallel do default(none), private(ato), shared(wfn0,wfn1,wfn2,matC,wfns0,wfns1,wfns2,obt,level,nmat)
        do ato=1,nmat ! 循环每一个原子轨道，计算分子轨道波函数
            if (abs(matC(obt,ato))<1e-4) cycle
            wfn0(:,obt) = wfn0(:,obt) + wfns0(:,ato)*matC(obt,ato)
            if (level==0) cycle
            wfn1(:,:,obt) = wfn1(:,:,obt) + wfns1(:,:,ato)*matC(obt,ato)
            if (level==1) cycle
            wfn2(:,:,:,obt) = wfn2(:,:,:,obt) + wfns2(:,:,:,ato)*matC(obt,ato)
        end do
        !$omp end parallel do

        den0(:,obt)=den0(:,obt) + wfn0(:,obt)*wfn0(:,obt) !电子密度是波函数的平方
        if (level==0) cycle
        do i=1,3
            den1(i,:,obt)=den1(i,:,obt)+2*wfn0(:,obt)*wfn1(i,:,obt)
        end do
        if (level==1) cycle
        do i=1,3
            do j=1,3
                den2(i,j,:,obt)=wfn1(i,:,obt)*wfn1(j,:,obt) + wfn0(:,obt)*wfn2(i,j,:,obt)
            end do
        end do
    end do
    dens0=sum(den0,dim=2)
    dens1=sum(den1,dim=3)
    dens2=sum(den2,dim=4)
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