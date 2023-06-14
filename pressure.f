c-----------------------------------------------------------------------
c
      subroutine pressure(fac)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ic,jc,kc
      real*8 fac,ux,uy,uz,vx,vy,vz,wx,wy,wz

      do k=1,nz
        do j=1,ny
          do i=1,nx
            d(i,j,k)=data1(i,j,k)+data5(i,j,k)+data9(i,j,k)
          enddo
        enddo
      enddo

      do k=1,nz
        do j=1,ny
          do i=1,nx
            ux=data1(i,j,k)
            uy=data2(i,j,k)
            uz=data3(i,j,k)
            vx=data4(i,j,k)
            vy=data5(i,j,k)
            vz=data6(i,j,k)
            wx=data7(i,j,k)
            wy=data8(i,j,k)
            wz=data9(i,j,k)
            o(i,j,k)=dn(i,j,k)/(fac*dt)-
     .               dx1*(u(i,j,k)*(d(i+1,j,k)+d(i,j,k))-
     .                    u(i-1,j,k)*(d(i,j,k)+d(i-1,j,k)))-
     .               dy1*(v(i,j,k)*(d(i,j+1,k)+d(i,j,k))-
     .                    v(i,j-1,k)*(d(i,j,k)+d(i,j-1,k)))-
     .               dz1*(w(i,j,k)*(d(i,j,k+1)+d(i,j,k))-
     .                    w(i,j,k-1)*(d(i,j,k)+d(i,j,k-1)))+Re1*
     .               (dx2*(d(i+1,j,k)-2.0d0*d(i,j,k)+d(i-1,j,k))+
     .                dy2*(d(i,j+1,k)-2.0d0*d(i,j,k)+d(i,j-1,k))+
     .                dz2*(d(i,j,k+1)-2.0d0*d(i,j,k)+d(i,j,k-1)))+
     .               2.0d0*(ux*vy-vx*uy+ux*wz-wx*uz+vy*wz-wy*vz)
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

      call pbc
      call poisson_fft

      do k=1,nz
        do j=1,ny
          do i=1,nx
            fccc(i,j,k)=p(i,j,k)
          enddo
        enddo
      enddo      

      return
      end


c-----------------------------------------------------------------------
