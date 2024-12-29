module march
    use iso_c_binding
    implicit none
contains
    subroutine grids2voxel(nx,ny,nz, ngrid, grids, values, voxels)
        integer(c_int), intent(in),value :: nx,ny,nz,ngrid
        real(c_double), intent(in) :: grids(3,ngrid), values(ngrid)
        real(c_double), intent(inout) :: voxels(4,nz,ny,nx)
        
        integer:: i,j,k,n
        n=1
        do i=1,nx
            do j=1,ny
                do k=1,nz
                    voxels(:3,k,j,i)=grids(:,n)
                    voxels( 4,k,j,i)=values(n)
                    n=n+1
                end do
            end do
        end do
        
    end subroutine grids2voxel

    ! subroutine voxel2verts(nx,ny,nz,voxel,minl,maxl,gt,verts) !输出的顶点数量是不固定的，这可咋整
    !     integer(c_int), intent(in),value :: nx,ny,nz,minl,maxl,gt
    !     real(c_double), intent(in):: voxel(4,nz,ny,nx)
    !     real(c_double), intent(inout):: verts
    
        
    ! end subroutine voxel2verts
end module march