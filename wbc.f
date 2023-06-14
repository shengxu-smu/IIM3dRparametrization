c-----------------------------------------------------------------------
c
      subroutine wbc
      include 'parameter.inc'
      include 'field.inc'

c  periodic
      if(lb_wbc.eq.0.and.lt_wbc.eq.0) then
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,0)=w(i,j,nz-1)
            w(i,j,nz)=w(i,j,1)
          enddo
        enddo
      endif

      if(ls_wbc.eq.0.and.ln_wbc.eq.0) then
        do k=0,nz
          do i=0,nx+1
            w(i,0,k)=w(i,ny-1,k)
            w(i,ny+1,k)=w(i,2,k)
          enddo
        enddo
      endif

      if(lw_wbc.eq.0.and.le_wbc.eq.0) then
        do k=0,nz
          do j=0,ny+1
            w(0,j,k)=w(nx-1,j,k)
            w(nx+1,j,k)=w(2,j,k)
          enddo
        enddo
      endif

c  dirichlet

      if(lb_wbc.eq.1) then
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,0)=-w(i,j,1)
          enddo
        enddo
      endif
      if(lt_wbc.eq.1) then
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,nz)=-w(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_wbc.eq.1) then
        do k=0,nz
          do i=0,nx+1
            w(i,0,k)=-w(i,2,k)
          enddo
        enddo
      endif
      if(ln_wbc.eq.1) then
        do k=0,nz
          do i=0,nx+1
            w(i,ny+1,k)=-w(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_wbc.eq.1) then
        do k=0,nz
          do j=0,ny+1
            w(0,j,k)=-w(2,j,k)
          enddo
        enddo
      endif
      if(le_wbc.eq.1) then
        do k=0,nz
          do j=0,ny+1
            w(nx+1,j,k)=-w(nx-1,j,k)
          enddo
        enddo
      endif

c  newmann

      if(lb_wbc.eq.2) then
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,0)=w(i,j,1)
          enddo
        enddo
      endif
      if(lt_wbc.eq.2) then
        do j=0,ny+1
          do i=0,nx+1
            w(i,j,nz)=w(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_wbc.eq.2) then
        do k=0,nz
          do i=0,nx+1
            w(i,0,k)=w(i,2,k)
          enddo
        enddo
      endif
      if(ln_wbc.eq.2) then
        do k=0,nz
          do i=0,nx+1
            w(i,ny+1,k)=w(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_wbc.eq.2) then
        do k=0,nz
          do j=0,ny+1
            w(0,j,k)=w(2,j,k)
          enddo
        enddo
      endif
      if(le_wbc.eq.2) then
        do k=0,nz
          do j=0,ny+1
            w(nx+1,j,k)=w(nx-1,j,k)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
