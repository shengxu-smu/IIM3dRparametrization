c-----------------------------------------------------------------------
c
      subroutine u_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx
            fecc(i,j,k)=u(i,j,k)
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny+1
          do i=0,nx
            fece(i,j,k)=0.5d0*(u(i,j,k)+u(i,j,k+1))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=0,nx
            feec(i,j,k)=0.5d0*(u(i,j,k)+u(i,j+1,k))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny+1
          do i=1,nx
            fccc(i,j,k)=0.5d0*(u(i,j,k)+u(i-1,j,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,niejc(m)
        i=iejc(1,n,m)
        j=iejc(2,n,m)
        k=kuiejc(n,m)
        fece(i,j,k)=fece(i,j,k)+ukiejc(n,m)
      enddo

      do n=0,niekc(m)
        i=iekc(1,n,m)
        k=iekc(2,n,m)
        j=juiekc(n,m)
        feec(i,j,k)=feec(i,j,k)+ujiekc(n,m)
      enddo

      do n=0,njckc(m)
        j=jckc(1,n,m)
        k=jckc(2,n,m)
        i=iujckc(n,m)
        fccc(i,j,k)=fccc(i,j,k)+uijckc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=0,nz
        do j=0,ny+1
          do i=1,nx
            fcce(i,j,k)=0.5d0*(fece(i,j,k)+fece(i-1,j,k))
          enddo
        enddo
      enddo

      do k=0,nz+1
        do j=0,ny
          do i=1,nx
            fcec(i,j,k)=0.5d0*(feec(i,j,k)+feec(i-1,j,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,njcke(m)
        j=jcke(1,n,m)
        k=jcke(2,n,m)
        i=iujcke(n,m)
        fcce(i,j,k)=fcce(i,j,k)+uijcke(n,m)
      enddo

      do n=0,njekc(m)
        j=jekc(1,n,m)
        k=jekc(2,n,m)
        i=iujekc(n,m)
        fcec(i,j,k)=fcec(i,j,k)+uijekc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=0,nz
        do j=0,ny
          do i=1,nx
            fcee(i,j,k)=0.5d0*(fcec(i,j,k)+fcec(i,j,k+1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicje(m)
        i=icje(1,n,m)
        j=icje(2,n,m)
        k=kuicje(n,m)
        fcee(i,j,k)=fcee(i,j,k)+ukicje(n,m)
      enddo

      ENDDO
      ENDIF

      return
      end


c-----------------------------------------------------------------------



