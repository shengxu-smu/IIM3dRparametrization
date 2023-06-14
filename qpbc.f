c-----------------------------------------------------------------------
c
      subroutine qpbc
      include 'parameter.inc'
      include 'field.inc'

      do j=1,ny
        do i=1,nx
          pb(i,j)=0.0d0
          pt(i,j)=0.0d0
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          pw(j,k)=0.0d0
          pe(j,k)=0.0d0
        enddo
        do i=1,nx
          ps(i,k)=0.0d0
          pn(i,k)=0.0d0
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
