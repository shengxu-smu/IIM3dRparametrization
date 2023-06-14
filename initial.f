c-----------------------------------------------------------------------
c
      subroutine initial
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      include 'fft.inc'
      real*8 tc,rx,ry,rz,rr,c,s,e
      real*8 rxxt,rxyt,rxzt,ryxt,ryyt,ryzt,rzxt,rzyt,rzzt

      call vrffti(ns1,wsave1)
      call vrffti(ns2,wsave2)
      call vrffti(2*ns1+2,wsave3)

      call vcosti(nx,wsavei)
      call vcosti(ny,wsavej)
      call vcosti(nz,wsavek)

      t=0.0d0
      t0=0.0d0

      open(unit=21,file='DAT/alfa10.dat',status='unknown')
      open(unit=22,file='DAT/alfa11.dat',status='unknown')
      open(unit=23,file='DAT/alfa20.dat',status='unknown')
      do m1=0,ns1
        alfa10(m1)=dalfa10*dble(m1)
        alfa11(m1)=dalfa11*dble(m1)+0.5d0*dalfa11
        write(21,100)alfa10(m1)
        write(22,100)alfa11(m1)
      enddo
      do m2=0,ns2
        alfa20(m2)=dalfa20*dble(m2)
        write(23,100)alfa20(m2)
      enddo
      close(21)
      close(22)
      close(23)

      call surface_parametrization
      call index_initialize
      call surface_initialize
      open(unit=24,file='DAT/xs0.dat',status='unknown')
      open(unit=25,file='DAT/ys0.dat',status='unknown')
      open(unit=26,file='DAT/zs0.dat',status='unknown')
      DO m=1,ms
      do m2=0,ns2
        write(24,100)(xs(m1,m2,m),m1=0,ns1)
        write(25,100)(ys(m1,m2,m),m1=0,ns1)
        write(26,100)(zs(m1,m2,m),m1=0,ns1)
      enddo
      ENDDO
      close(24)
      close(25)
      close(26)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx
            u(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=0,nx+1
            v(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            p(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            d(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      call velocity_reset

100   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
