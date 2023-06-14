c-----------------------------------------------------------------------
c
      subroutine qpressure
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ic,jc,kc
      real*8 fac,ux,uy,uz,vx,vy,vz,wx,wy,wz

      do k=1,nz
        do j=1,ny
          do i=1,nx
            o(i,j,k)=0.0d0
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,njckc(m)
          jc=jckc(1,n,m)
          kc=jckc(2,n,m)
          ic=jckc(4,n,m)
          o(ic,jc,kc)=o(ic,jc,kc)-pdxx(1,n,m)
          o(ic+1,jc,kc)=o(ic+1,jc,kc)-pdxx(2,n,m)
        enddo
        do n=0,nickc(m)
          ic=ickc(1,n,m)
          kc=ickc(2,n,m)
          jc=ickc(4,n,m)
          o(ic,jc,kc)=o(ic,jc,kc)-pdyy(1,n,m)
          o(ic,jc+1,kc)=o(ic,jc+1,kc)-pdyy(2,n,m)
        enddo
        do n=0,nicjc(m)
          ic=icjc(1,n,m)
          jc=icjc(2,n,m)
          kc=icjc(4,n,m)
          o(ic,jc,kc)=o(ic,jc,kc)-pdzz(1,n,m)
          o(ic,jc,kc+1)=o(ic,jc,kc+1)-pdzz(2,n,m)
        enddo
        ENDDO
      endif

      call qpbc
      call poisson_fft

      return
      end


c-----------------------------------------------------------------------
