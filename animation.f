c-----------------------------------------------------------------------
c
      subroutine animation
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      call vorticity

      DO m=1,ms
      do m2=0,ns2
        write(63,200)(xs(m1,m2,m),m1=0,ns1)
        write(64,200)(ys(m1,m2,m),m1=0,ns1)
        write(65,200)(zs(m1,m2,m),m1=0,ns1)
      enddo
      ENDDO

      do k=1,nz
        do j=1,ny
          write(66,200)(o(i,j,k),i=1,nx)
        enddo
      enddo

200   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
