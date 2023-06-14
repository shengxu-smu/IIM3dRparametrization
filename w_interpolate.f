c-----------------------------------------------------------------------
c
      subroutine w_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      do k=0,nz
        do j=0,ny+1
          do i=0,nx+1
            fcce(i,j,k)=w(i,j,k)
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny
          do i=0,nx+1
            fcee(i,j,k)=0.5d0*(w(i,j,k)+w(i,j+1,k))
          enddo
        enddo
      enddo

      do k=0,nz
        do j=0,ny+1
          do i=0,nx
            fece(i,j,k)=0.5d0*(w(i,j,k)+w(i+1,j,k))
          enddo
        enddo
      enddo

      do k=1,nz
        do j=0,ny+1
          do i=0,nx+1
            fccc(i,j,k)=0.5d0*(w(i,j,k)+w(i,j,k-1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicke(m)
        i=icke(1,n,m)
        k=icke(2,n,m)
        j=jwicke(n,m)
        fcee(i,j,k)=fcee(i,j,k)+wjicke(n,m)
      enddo

      do n=0,njcke(m)
        j=jcke(1,n,m)
        k=jcke(2,n,m)
        i=iwjcke(n,m)
        fece(i,j,k)=fece(i,j,k)+wijcke(n,m)
      enddo

      do n=0,nicjc(m)
        i=icjc(1,n,m)
        j=icjc(2,n,m)
        k=kwicjc(n,m)
        fccc(i,j,k)=fccc(i,j,k)+wkicjc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=1,nz
        do j=0,ny
          do i=0,nx+1
            fcec(i,j,k)=0.5d0*(fcee(i,j,k)+fcee(i,j,k-1))
          enddo
        enddo
      enddo

      do k=1,nz
        do j=0,ny+1
          do i=0,nx
            fecc(i,j,k)=0.5d0*(fece(i,j,k)+fece(i,j,k-1))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN
      DO m=1,ms

      do n=0,nicje(m)
        i=icje(1,n,m)
        j=icje(2,n,m)
        k=kwicje(n,m)
        fcec(i,j,k)=fcec(i,j,k)+wkicje(n,m)
      enddo

      do n=0,niejc(m)
        i=iejc(1,n,m)
        j=iejc(2,n,m)
        k=kwiejc(n,m)
        fecc(i,j,k)=fecc(i,j,k)+wkiejc(n,m)
      enddo

      ENDDO
      ENDIF

      do k=1,nz
        do j=0,ny
          do i=0,nx
            feec(i,j,k)=0.5d0*(fecc(i,j,k)+fecc(i,j+1,k))
          enddo
        enddo
      enddo

      IF(isingular.eq.1) THEN     
      DO m=1,ms

      do n=0,niekc(m)
        i=iekc(1,n,m)
        k=iekc(2,n,m)
        j=jwiekc(n,m)
        feec(i,j,k)=feec(i,j,k)+wjiekc(n,m)
      enddo

      ENDDO
      ENDIF

      return
      end


c-----------------------------------------------------------------------



