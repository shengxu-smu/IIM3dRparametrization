c-----------------------------------------------------------------------
c
      subroutine ubc
      include 'parameter.inc'
      include 'field.inc'

c  periodic
      if(lb_ubc.eq.0.and.lt_ubc.eq.0) then
        do j=0,ny+1
          do i=0,nx
            u(i,j,0)=u(i,j,nz-1)
            u(i,j,nz+1)=u(i,j,2)
          enddo
        enddo
      endif

      if(ls_ubc.eq.0.and.ln_ubc.eq.0) then
        do k=0,nz+1
          do i=0,nx
            u(i,0,k)=u(i,ny-1,k)
            u(i,ny+1,k)=u(i,2,k)
          enddo
        enddo
      endif

      if(lw_ubc.eq.0.and.le_ubc.eq.0) then
        do k=0,nz+1
          do j=0,ny+1
            u(0,j,k)=u(nx-1,j,k)
            u(nx,j,k)=u(1,j,k)
          enddo
        enddo
      endif

c  dirichlet

      if(lb_ubc.eq.1) then
        do j=0,ny+1
          do i=0,nx
            u(i,j,0)=0.0d0-u(i,j,2)
          enddo
        enddo
      endif
      if(lt_ubc.eq.1) then
        do j=0,ny+1
          do i=0,nx
            u(i,j,nz+1)=0.0d0-u(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_ubc.eq.1) then
        do k=0,nz+1
          do i=0,nx
            u(i,0,k)=0.0d0-u(i,2,k)
          enddo
        enddo
      endif
      if(ln_ubc.eq.1) then
        do k=0,nz+1
          do i=0,nx
            u(i,ny+1,k)=0.0d0-u(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_ubc.eq.1) then
        do k=0,nz+1
          do j=0,ny+1
            u(0,j,k)=0.0d0-u(1,j,k)
          enddo
        enddo
      endif
      if(le_ubc.eq.1) then
        do k=0,nz+1
          do j=0,ny+1
            u(nx,j,k)=0.0d0-u(nx-1,j,k)
          enddo
        enddo
      endif

c  newmann

      if(lb_ubc.eq.2) then
        do j=0,ny+1
          do i=0,nx
            u(i,j,0)=u(i,j,2)
          enddo
        enddo
      endif
      if(lt_ubc.eq.2) then
        do j=0,ny+1
          do i=0,nx
            u(i,j,nz+1)=u(i,j,nz-1)
          enddo
        enddo
      endif

      if(ls_ubc.eq.2) then
        do k=0,nz+1
          do i=0,nx
            u(i,0,k)=u(i,2,k)
          enddo
        enddo
      endif
      if(ln_ubc.eq.2) then
        do k=0,nz+1
          do i=0,nx
            u(i,ny+1,k)=u(i,ny-1,k)
          enddo
        enddo
      endif

      if(lw_ubc.eq.2) then
        do k=0,nz+1
          do j=0,ny+1
            u(0,j,k)=u(1,j,k)
          enddo
        enddo
      endif
      if(le_ubc.eq.2) then
        do k=0,nz+1
          do j=0,ny+1
            u(nx,j,k)=u(nx-1,j,k)
          enddo
        enddo
      endif

      return
      end


c-----------------------------------------------------------------------
