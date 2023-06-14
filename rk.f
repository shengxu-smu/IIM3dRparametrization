c-----------------------------------------------------------------------
c
      subroutine rk
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer krk,modify
      real*8 fac

      modify=0

c  krk=1:
      krk=1
      fac=0.5d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      call old_save
      call pressure(fac)
      call rhs

      if(modify.eq.1) then
        call surface_pressure
        call gmres_pressure
        call velocity_correct
      endif

      if(imove.eq.1) call surface_move

!$OMP PARALLEL PRIVATE(I,J,K)

!$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
            um(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
            vm(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
            wm(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      call velocity_reset
      call ubc
      call vbc
      call wbc

c  krk=2:
      krk=2
      fac=0.5d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      call pressure(fac)
      call rhs

      if(modify.eq.1) then
        call surface_pressure
        call gmres_pressure
        call velocity_correct
      endif

!$OMP PARALLEL PRIVATE(I,J,K)
!$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
            um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
            vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
            wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      call velocity_reset
      call ubc
      call vbc
      call wbc

c  krk=3:
      krk=3
      fac=1.0d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      call pressure(fac)
      call rhs

      if(modify.eq.1) then
        call surface_pressure
        call gmres_pressure
        call velocity_correct
      endif

      if(imove.eq.1) call surface_move

!$OMP PARALLEL PRIVATE(I,J,K)
!$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k)=un(i,j,k)+fac*dt*data1(i,j,k)
            v(i,j,k)=vn(i,j,k)+fac*dt*data4(i,j,k)
            w(i,j,k)=wn(i,j,k)+fac*dt*data7(i,j,k)
            um(i,j,k)=um(i,j,k)+(dt/3.0d0)*data1(i,j,k)
            vm(i,j,k)=vm(i,j,k)+(dt/3.0d0)*data4(i,j,k)
            wm(i,j,k)=wm(i,j,k)+(dt/3.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      call velocity_reset
      call ubc
      call vbc
      call wbc

c  krk=4:
      krk=4
      fac=1.0d0

      if(isingular.eq.1) then
        call singular_call
      else
        call u_interpolate
        call u_strain
        call v_interpolate
        call v_strain
        call w_interpolate
        call w_strain
      endif

      call pressure(fac)
      call rhs

      if(modify.eq.1) then
        call surface_pressure
        call gmres_pressure
        call velocity_correct
      endif

!$OMP PARALLEL PRIVATE(I,J,K)
!$OMP DO
      do k=1,nz
        do j=1,ny
          do i=1,nx
            u(i,j,k)=um(i,j,k)+(dt/6.0d0)*data1(i,j,k)
            v(i,j,k)=vm(i,j,k)+(dt/6.0d0)*data4(i,j,k)
            w(i,j,k)=wm(i,j,k)+(dt/6.0d0)*data7(i,j,k)
          enddo
        enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL

      call velocity_reset
      call ubc
      call vbc
      call wbc

      return
      end


c-----------------------------------------------------------------------
