c-----------------------------------------------------------------------
c
      subroutine rhs
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc
      real*8 conv,grad,visc

      do k=1,nz
        do j=1,ny
          do i=1,nx-1
            conv=-dx1*(data1(i+1,j,k)*data1(i+1,j,k)-
     .                 data1(i,j,k)*data1(i,j,k))
     .           -dy1*(data2(i,j,k)*data5(i,j,k)-
     .                 data2(i,j-1,k)*data5(i,j-1,k))
     .           -dz1*(data3(i,j,k)*data8(i,j,k)-
     .                 data3(i,j,k-1)*data8(i,j,k-1))
            grad=-dx1*(p(i+1,j,k)-p(i,j,k))
            visc=dx2*(u(i+1,j,k)-2.0d0*u(i,j,k)+u(i-1,j,k))+
     .           dy2*(u(i,j+1,k)-2.0d0*u(i,j,k)+u(i,j-1,k))+
     .           dz2*(u(i,j,k+1)-2.0d0*u(i,j,k)+u(i,j,k-1))
            data1(i,j,k)=conv+grad+Re1*visc
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,njckc(m)
          jc=jckc(1,n,m)
          kc=jckc(2,n,m)
          ie=jckc(3,n,m)
          ic=jckc(4,n,m)
          data1(ic,jc,kc)=data1(ic,jc,kc)-uudx(n,m)-pdx(n,m)
          data1(ie,jc,kc)=data1(ie,jc,kc)+Re1*udxx(1,n,m)
          data1(ie+1,jc,kc)=data1(ie+1,jc,kc)+Re1*udxx(2,n,m)
        enddo
        do n=0,niekc(m)
          ie=iekc(1,n,m)
          kc=iekc(2,n,m)
          je=iekc(3,n,m)
          jc=iekc(4,n,m)
          data1(ie,je+1,kc)=data1(ie,je+1,kc)-uvdy(n,m)
          data1(ie,jc,kc)=data1(ie,jc,kc)+Re1*udyy(1,n,m)
          data1(ie,jc+1,kc)=data1(ie,jc+1,kc)+Re1*udyy(2,n,m)
        enddo
        do n=0,niejc(m)
          ie=iejc(1,n,m)
          jc=iejc(2,n,m)
          ke=iejc(3,n,m)
          kc=iejc(4,n,m)
          data1(ie,jc,ke+1)=data1(ie,jc,ke+1)-uwdz(n,m)
          data1(ie,jc,kc)=data1(ie,jc,kc)+Re1*udzz(1,n,m)
          data1(ie,jc,kc+1)=data1(ie,jc,kc+1)+Re1*udzz(2,n,m)
        enddo
        ENDDO
      endif

      do k=1,nz
        do j=1,ny-1
          do i=1,nx
            conv=-dx1*(data5(i,j,k)*data2(i,j,k)-
     .                 data5(i-1,j,k)*data2(i-1,j,k))
     .           -dy1*(data4(i,j+1,k)*data4(i,j+1,k)-
     .                 data4(i,j,k)*data4(i,j,k))
     .           -dz1*(data6(i,j,k)*data9(i,j,k)-
     .                 data6(i,j,k-1)*data9(i,j,k-1))
            grad=-dy1*(p(i,j+1,k)-p(i,j,k))
            visc=dx2*(v(i+1,j,k)-2.0d0*v(i,j,k)+v(i-1,j,k))+
     .           dy2*(v(i,j+1,k)-2.0d0*v(i,j,k)+v(i,j-1,k))+
     .           dz2*(v(i,j,k+1)-2.0d0*v(i,j,k)+v(i,j,k-1))
            data4(i,j,k)=conv+grad+Re1*visc
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,nickc(m)
          ic=ickc(1,n,m)
          kc=ickc(2,n,m)
          je=ickc(3,n,m)
          jc=ickc(4,n,m)
          data4(ic,jc,kc)=data4(ic,jc,kc)-vvdy(n,m)-pdy(n,m)
          data4(ic,je,kc)=data4(ic,je,kc)+Re1*vdyy(1,n,m)
          data4(ic,je+1,kc)=data4(ic,je+1,kc)+Re1*vdyy(2,n,m)
        enddo
        do n=0,njekc(m)
          je=jekc(1,n,m)
          kc=jekc(2,n,m)
          ie=jekc(3,n,m)
          ic=jekc(4,n,m)
          data4(ie+1,je,kc)=data4(ie+1,je,kc)-vudx(n,m)
          data4(ic,je,kc)=data4(ic,je,kc)+Re1*vdxx(1,n,m)
          data4(ic+1,je,kc)=data4(ic+1,je,kc)+Re1*vdxx(2,n,m)
        enddo
        do n=0,nicje(m)
          ic=icje(1,n,m)
          je=icje(2,n,m)
          ke=icje(3,n,m)
          kc=icje(4,n,m)
          data4(ic,je,ke+1)=data4(ic,je,ke+1)-vwdz(n,m)
          data4(ic,je,kc)=data4(ic,je,kc)+Re1*vdzz(1,n,m)
          data4(ic,je,kc+1)=data4(ic,je,kc+1)+Re1*vdzz(2,n,m)
        enddo
        ENDDO
      endif

      do k=1,nz-1
        do j=1,ny
          do i=1,nx
            conv=-dx1*(data8(i,j,k)*data3(i,j,k)-
     .                 data8(i-1,j,k)*data3(i-1,j,k))
     .           -dy1*(data9(i,j,k)*data6(i,j,k)-
     .                 data9(i,j-1,k)*data6(i,j-1,k))
     .           -dz1*(data7(i,j,k+1)*data7(i,j,k+1)-
     .                 data7(i,j,k)*data7(i,j,k))
            grad=-dz1*(p(i,j,k+1)-p(i,j,k))
            visc=dx2*(w(i+1,j,k)-2.0d0*w(i,j,k)+w(i-1,j,k))+
     .           dy2*(w(i,j+1,k)-2.0d0*w(i,j,k)+w(i,j-1,k))+
     .           dz2*(w(i,j,k+1)-2.0d0*w(i,j,k)+w(i,j,k-1))
            data7(i,j,k)=conv+grad+Re1*visc
          enddo
        enddo
      enddo

      if(isingular.eq.1) then
        DO m=1,ms
        do n=0,nicjc(m)
          ic=icjc(1,n,m)
          jc=icjc(2,n,m)
          ke=icjc(3,n,m)
          kc=icjc(4,n,m)
          data7(ic,jc,kc)=data7(ic,jc,kc)-wwdz(n,m)-pdz(n,m)
          data7(ic,jc,ke)=data7(ic,jc,ke)+Re1*wdzz(1,n,m)
          data7(ic,jc,ke+1)=data7(ic,jc,ke+1)+Re1*wdzz(2,n,m)
        enddo
        do n=0,njcke(m)
          jc=jcke(1,n,m)
          ke=jcke(2,n,m)
          ie=jcke(3,n,m)
          ic=jcke(4,n,m)
          data7(ie+1,jc,ke)=data7(ie+1,jc,ke)-wudx(n,m)
          data7(ic,jc,ke)=data7(ic,jc,ke)+Re1*wdxx(1,n,m)
          data7(ic+1,jc,ke)=data7(ic+1,jc,ke)+Re1*wdxx(2,n,m)
        enddo
        do n=0,nicke(m)
          ic=icke(1,n,m)
          ke=icke(2,n,m)
          je=icke(3,n,m)
          jc=icke(4,n,m)
          data7(ic,je+1,ke)=data7(ic,je+1,ke)-wvdy(n,m)
          data7(ic,jc,ke)=data7(ic,jc,ke)+Re1*wdyy(1,n,m)
          data7(ic,jc+1,ke)=data7(ic,jc+1,ke)+Re1*wdyy(2,n,m)
        enddo
        ENDDO
      endif

      return
      end


c-----------------------------------------------------------------------
