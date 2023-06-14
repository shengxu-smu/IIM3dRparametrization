c-----------------------------------------------------------------------
c
      subroutine w_strain
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke

      do k=1,nz
        do j=1,ny
          do i=1,nx
            data7(i,j,k)=dx1*(fecc(i,j,k)-fecc(i-1,j,k))
            data8(i,j,k)=dy1*(fcec(i,j,k)-fcec(i,j-1,k))
            data9(i,j,k)=dz1*(fcce(i,j,k)-fcce(i,j,k-1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicjc(m)
        i=icjc(1,n,m)
        j=icjc(2,n,m)
        ke=icjc(3,n,m)
        k=ke+1
        data9(i,j,k)=data9(i,j,k)+wdz(n,m)
      enddo

      do n=0,nickc(m)
        i=ickc(1,n,m)
        k=ickc(2,n,m)
        je=ickc(3,n,m)
        j=je+1
        data8(i,j,k)=data8(i,j,k)+wdy(n,m)
      enddo

      do n=0,njckc(m)
        j=jckc(1,n,m)
        k=jckc(2,n,m)
        ie=jckc(3,n,m)
        i=ie+1
        data7(i,j,k)=data7(i,j,k)+wdx(n,m)
      enddo

      ENDDO
      ENDIF

      return
      end


c-----------------------------------------------------------------------



