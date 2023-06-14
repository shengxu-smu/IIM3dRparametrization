c-----------------------------------------------------------------------
c
      subroutine old_save
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'

!$OMP PARALLEL PRIVATE(I,J,K)
!$OMP DO
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx
            un(i,j,k)=u(i,j,k)
            um(i,j,k)=u(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO

!$OMP DO 
      do k=0,nz+1
        do j=0,ny
          do i=0,nx+1
            vn(i,j,k)=v(i,j,k)
            vm(i,j,k)=v(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO

!$OMP DO
      do k=0,nz
        do j=0,ny+1
          do i=0,nx+1
            wn(i,j,k)=w(i,j,k)
            wm(i,j,k)=w(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO

!$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1,nx
            dn(i,j,k)=data1(i,j,k)+data5(i,j,k)+data9(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      return
      end


c-----------------------------------------------------------------------

