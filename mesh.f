c-----------------------------------------------------------------------
c
      subroutine mesh
      include 'parameter.inc'
      include 'field.inc'

      dx=xl/dble(nx)
      hdx=0.5d0*dx
      dx1=1.0d0/dx
      dx2=dx1*dx1

      dy=yl/dble(ny)
      hdy=0.5d0*dy
      dy1=1.0d0/dy
      dy2=dy1*dy1

      dz=zl/dble(nz)
      hdz=0.5d0*dz
      dz1=1.0d0/dz
      dz2=dz1*dz1

      dalfa10=2.0d0*sl1/dble(ns1)
      dalfa11=sl1/dble(ns1+1)
      dalfa20=sl2/dble(ns2)

      betay=dx*dx*dy2
      betaz=dx*dx*dz2

      do i=-2,2*nx
        x(i)=x0+hdx*dble(i)
      enddo
      do j=-2,2*ny
        y(j)=y0+hdy*dble(j)
      enddo
      do k=-2,2*nz
        z(k)=z0+hdz*dble(k)
      enddo

c staggered MAC grid

      xc(0)=x(-2)
      do i=0,nx
        xc(i+1)=x(2*i)
        xe(i)=xc(i)+hdx
      enddo
      yc(0)=y(-2)
      do j=0,ny
        yc(j+1)=y(2*j)
        ye(j)=yc(j)+hdy
      enddo
      zc(0)=z(-2)
      do k=0,nz
        zc(k+1)=z(2*k)
        ze(k)=zc(k)+hdz
      enddo

      return
      end


c-----------------------------------------------------------------------
