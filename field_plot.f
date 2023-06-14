c-----------------------------------------------------------------------
c
      subroutine field_plot
      include 'parameter.inc'
      include 'field.inc'
      integer i0,j0,k0

      open(unit=81,file='DAT/xc.dat',status='unknown')
      open(unit=82,file='DAT/yc.dat',status='unknown')
      open(unit=83,file='DAT/zc.dat',status='unknown')
      do i=1,nx
        write(81,100)xc(i)
      enddo
      do j=1,ny
        write(82,100)yc(j)
      enddo
      do k=1,nz
        write(83,100)zc(k)
      enddo
      close(81)
      close(82)
      close(83)

      call vorticity

      i0=1-x0/dx
      j0=1-y0/dy
      k0=1-z0/dz
      open(unit=84,file='DAT/xu.dat',status='unknown')
      open(unit=85,file='DAT/yv.dat',status='unknown')
      open(unit=86,file='DAT/zw.dat',status='unknown')
      open(unit=87,file='DAT/xp.dat',status='unknown')
      open(unit=88,file='DAT/yp.dat',status='unknown')
      open(unit=89,file='DAT/zp.dat',status='unknown')
      do i=1,nx
        write(84,100)data1(i,j0,k0)
        write(87,100)p(i,j0,k0)
      enddo
      do j=1,ny
        write(85,100)data4(i0,j,k0)
        write(88,100)p(i0,j,k0)
      enddo
      do k=1,nz
        write(86,100)data7(i0,j0,k)
        write(89,100)p(i0,j0,k)
      enddo
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)

100   format(1x,e16.6e4)

      open(unit=91,file='DAT/u.dat',status='unknown')
      open(unit=92,file='DAT/v.dat',status='unknown')
      open(unit=93,file='DAT/w.dat',status='unknown')
      open(unit=94,file='DAT/o1.dat',status='unknown')
      open(unit=95,file='DAT/o2.dat',status='unknown')
      open(unit=96,file='DAT/o3.dat',status='unknown')
      open(unit=97,file='DAT/p.dat',status='unknown')
      open(unit=98,file='DAT/d.dat',status='unknown')
      open(unit=99,file='DAT/q.dat',status='unknown')
      do k=1,nz
        do j=1,ny
          write(91,200)(data1(i,j,k),i=1,nx)
          write(92,200)(data4(i,j,k),i=1,nx)
          write(93,200)(data7(i,j,k),i=1,nx)
          write(94,200)(fecc(i,j,k),i=1,nx)
          write(95,200)(fcec(i,j,k),i=1,nx)
          write(96,200)(fcce(i,j,k),i=1,nx)
          write(97,200)(p(i,j,k),i=1,nx)
          write(98,200)(d(i,j,k),i=1,nx)
          write(99,200)(o(i,j,k),i=1,nx)
        enddo
      enddo
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(96)
      close(97)
      close(98)
      close(99)

200   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
