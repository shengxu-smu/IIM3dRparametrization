c-----------------------------------------------------------------------
c
      subroutine v_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do k=0,nz+1
        do j=0,ny
          do i=0,nx+1
            fcec(i,j,k)=v(i,j,k)
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny
          do i=0,nx+1
            fcee(i,j,k)=0.5d0*(v(i,j,k)+v(i,j,k+1))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=0,nx
            feec(i,j,k)=0.5d0*(v(i,j,k)+v(i+1,j,k))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=1,ny
          do i=0,nx+1
            fccc(i,j,k)=0.5d0*(v(i,j,k)+v(i,j-1,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicje(m)
        i=icje(1,n,m)
        j=icje(2,n,m)
        k=kvicje(n,m)
        fcee(i,j,k)=fcee(i,j,k)+vkicje(n,m)
      enddo

      do n=0,njekc(m)
        j=jekc(1,n,m)
        k=jekc(2,n,m)
        i=ivjekc(n,m)
        feec(i,j,k)=feec(i,j,k)+vijekc(n,m)
      enddo

      do n=0,nickc(m)
        i=ickc(1,n,m)
        k=ickc(2,n,m)
        j=jvickc(n,m)
        fccc(i,j,k)=fccc(i,j,k)+vjickc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=0,nz
        do j=1,ny
          do i=0,nx+1
            fcce(i,j,k)=0.5d0*(fcee(i,j,k)+fcee(i,j-1,k))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=1,ny
          do i=0,nx
            fecc(i,j,k)=0.5d0*(feec(i,j,k)+feec(i,j-1,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicke(m)
        i=icke(1,n,m)
        k=icke(2,n,m)
        j=jvicke(n,m)
        fcce(i,j,k)=fcce(i,j,k)+vjicke(n,m)
      enddo

      do n=0,niekc(m)
        i=iekc(1,n,m)
        k=iekc(2,n,m)
        j=jviekc(n,m)
        fecc(i,j,k)=fecc(i,j,k)+vjiekc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=0,nz
        do j=1,ny
          do i=0,nx
            fece(i,j,k)=0.5d0*(fcce(i,j,k)+fcce(i+1,j,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,njcke(m)
        j=jcke(1,n,m)
        k=jcke(2,n,m)
        i=ivjcke(n,m)
        fece(i,j,k)=fece(i,j,k)+vijcke(n,m)
      enddo

      ENDDO
      ENDIF

      return
      end


c-----------------------------------------------------------------------



