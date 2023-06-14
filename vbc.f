c-----------------------------------------------------------------------
c
      subroutine vbc
      include 'parameter.inc'
      include 'field.inc'

c  periodic
      if(lb_vbc.eq.0.and.lt_vbc.eq.0) then
        do j=0,ny
          do i=0,nx+1
            v(i,j,0)=v(i,j,nz-1)
            v(i,j,nz+1)=v(i,j,2)
          enddo
        enddo
      endif

      if(ls_vbc.eq.0.and.ln_vbc.eq.0) then
        do k=0,nz+1
          do i=0,nx+1
            u(i,0,k)=u(i,ny-1,k)
            u(i,ny,k)=u(i,1,k)
          enddo
        enddo
      endif

      if(lw_vbc.eq.0.and.le_vbc.eq.0) then
        do k=0,nz+1
          do j=0,ny
            v(0,j,k)=v(nx-1,j,k)
            v(nx+1,j,k)=v(2,j,k)
          enddo
        enddo
      endif

c  dirichlet

      if(lb_vbc.eq.1) then
        do j=0,ny
          do i=0,nx+1
            v(i,j,0)=-v(i,j,2)
          enddo
        enddo
      endif
      if(lt_vbc.eq.1) then
        do j=0,ny
          do i=0,nx+1
            v(i,j,nz+1)=-v(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_vbc.eq.1) then
        do k=0,nz+1
          do i=0,nx+1
            v(i,0,k)=-v(i,1,k)
          enddo
        enddo
      endif
      if(ln_vbc.eq.1) then
        do k=0,nz+1
          do i=0,nx+1
            v(i,ny,k)=-v(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_vbc.eq.1) then
        do k=0,nz+1
          do j=0,ny
            v(0,j,k)=-v(2,j,k)
          enddo
        enddo
      endif
      if(le_vbc.eq.1) then
        do k=0,nz+1
          do j=0,ny
            v(nx+1,j,k)=-v(nx-1,j,k)
          enddo
        enddo
      endif

c  newmann

      if(lb_vbc.eq.2) then
        do j=0,ny
          do i=0,nx+1
            v(i,j,0)=v(i,j,2)
          enddo
        enddo
      endif
      if(lt_vbc.eq.2) then
        do j=0,ny
          do i=0,nx+1
            v(i,j,nz+1)=v(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_vbc.eq.2) then
        do k=0,nz+1
          do i=0,nx+1
            v(i,0,k)=v(i,1,k)
          enddo
        enddo
      endif
      if(ln_vbc.eq.2) then
        do k=0,nz+1
          do i=0,nx+1
            v(i,ny,k)=v(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_vbc.eq.2) then
        do k=0,nz+1
          do j=0,ny
            v(0,j,k)=v(2,j,k)
          enddo
        enddo
      endif
      if(le_vbc.eq.2) then
        do k=0,nz+1
          do j=0,ny
            v(nx+1,j,k)=v(nx-1,j,k)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
