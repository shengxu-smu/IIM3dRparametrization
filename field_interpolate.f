c-----------------------------------------------------------------------
c
      subroutine field_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do k=0,nz
        do j=0,ny
          do i=1,nx
            data1(i,j,k)=0.5d0*(u(i,j,k)+u(i-1,j,k))
          enddo
        enddo
      enddo
      do k=0,nz
        do j=0,ny
          do i=0,nx
            data2(i,j,k)=0.5d0*(u(i,j,k)+u(i,j+1,k))
            data3(i,j,k)=0.5d0*(u(i,j,k)+u(i,j,k+1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,niejc(m)
        i=iejc(1,n,m)
        j=iejc(2,n,m)
        k=kuiejc(n,m)
        data3(i,j,k)=data3(i,j,k)+ukiejc(n,m)
      enddo

      do n=0,niekc(m)
        i=iekc(1,n,m)
        k=iekc(2,n,m)
        j=juiekc(n,m)
        data2(i,j,k)=data2(i,j,k)+ujiekc(n,m)
      enddo

      do n=0,njckc(m)
        j=jckc(1,n,m)
        k=jckc(2,n,m)
        i=iujckc(n,m)
        data1(i,j,k)=data1(i,j,k)+uijckc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=0,nz
        do j=1,ny
          do i=0,nx
            data4(i,j,k)=0.5d0*(v(i,j,k)+v(i,j-1,k))
          enddo
        enddo
      enddo
      do k=0,nz
        do j=0,ny
          do i=0,nx
            data5(i,j,k)=0.5d0*(v(i,j,k)+v(i+1,j,k))
            data6(i,j,k)=0.5d0*(v(i,j,k)+v(i,j,k+1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicje(m)
        i=icje(1,n,m)
        j=icje(2,n,m)
        k=kvicje(n,m)
        data6(i,j,k)=data6(i,j,k)+vkicje(n,m)
      enddo

      do n=0,njekc(m)
        j=jekc(1,n,m)
        k=jekc(2,n,m)
        i=ivjekc(n,m)
        data5(i,j,k)=data5(i,j,k)+vijekc(n,m)
      enddo

      do n=0,nickc(m)
        i=ickc(1,n,m)
        k=ickc(2,n,m)
        j=jvickc(n,m)
        data4(i,j,k)=data4(i,j,k)+vjickc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=1,nz
        do j=0,ny
          do i=0,nx
            data7(i,j,k)=0.5d0*(w(i,j,k)+w(i,j,k-1))
          enddo
        enddo
      enddo
      do k=0,nz
        do j=0,ny
          do i=0,nx
            data8(i,j,k)=0.5d0*(w(i,j,k)+w(i+1,j,k))
            data9(i,j,k)=0.5d0*(w(i,j,k)+w(i,j+1,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicke(m)
        i=icke(1,n,m)
        k=icke(2,n,m)
        j=jwicke(n,m)
        data9(i,j,k)=data9(i,j,k)+wjicke(n,m)
      enddo

      do n=0,njcke(m)
        j=jcke(1,n,m)
        k=jcke(2,n,m)
        i=iwjcke(n,m)
        data8(i,j,k)=data8(i,j,k)+wijcke(n,m)
      enddo
      do n=0,nicjc(m)
        i=icjc(1,n,m)
        j=icjc(2,n,m)
        k=kwicjc(n,m)
        data7(i,j,k)=data7(i,j,k)+wkicjc(n,m)
      enddo

      ENDDO
      ENDIF

      return
      end


c-----------------------------------------------------------------------



