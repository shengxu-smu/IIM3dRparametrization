c-----------------------------------------------------------------------
c
      subroutine velocity_correct
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc
      real*8 conv,grad,visc

      call qjc_rhs
      call qjc_first
      call qjc_second
      call qeuler_link
      call qjc_pressure
      call qcorrection_difference
      call qpressure

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            grad=-dx1*(p(i+1,j,k)-p(i,j,k))
            data1(i,j,k)=data1(i,j,k)+grad
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,njckc(m)
          jc=jckc(1,n,m)
          kc=jckc(2,n,m)
          ic=jckc(4,n,m)
          data1(ic,jc,kc)=data1(ic,jc,kc)-pdx(n,m)
        enddo
        ENDDO
      endif

      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            grad=-dy1*(p(i,j+1,k)-p(i,j,k))
            data4(i,j,k)=data4(i,j,k)+grad
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,nickc(m)
          ic=ickc(1,n,m)
          kc=ickc(2,n,m)
          jc=ickc(4,n,m)
          data4(ic,jc,kc)=data4(ic,jc,kc)-pdy(n,m)
        enddo
        ENDDO
      endif

      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            grad=-dz1*(p(i,j,k+1)-p(i,j,k))
            data7(i,j,k)=data7(i,j,k)+grad
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,nicjc(m)
          ic=icjc(1,n,m)
          jc=icjc(2,n,m)
          kc=icjc(4,n,m)
          data7(ic,jc,kc)=data7(ic,jc,kc)-pdz(n,m)
        enddo
        ENDDO
      endif

      do k=1,nz
        do j=1,ny
          do i=1,nx
            p(i,j,k)=fccc(i,j,k)+p(i,j,k)
          enddo
        enddo
      enddo

      DO m=1,ms
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,4,m)=fr(m1,m2,4,m)+fp(m1,m2,m)
        enddo
      enddo
      ENDDO

      return
      end


c-----------------------------------------------------------------------
