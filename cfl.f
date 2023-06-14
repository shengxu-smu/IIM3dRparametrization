c-----------------------------------------------------------------------
c
      subroutine cfl
      include 'parameter.inc'
      include 'field.inc'
      real*8 umax,vmax,wmax,dtc,dtv

      if(icfl.eq.0) then
        dt=dtfix
        t=t+dt
        return
      endif

      umax=0.0d0
      vmax=0.0d0
      wmax=0.0d0
!$OMP PARALLEL DO PRIVATE(I,J,K) REDUCTION(max:umax,vmax,wmax)
      do k=1,nz
        do j=1,ny
          do i=1,nx
            umax=max(umax,abs(u(i,j,k)))
            vmax=max(vmax,abs(v(i,j,k)))
            wmax=max(wmax,abs(w(i,j,k)))
          enddo
        enddo
      enddo
!OMP END PARALLEL DO

      if(umax.eq.0.0d0.and.vmax.eq.0.0d0.and.wmax.eq.0.0d0) then
        dtc=100.0d0
      else
        dtc=cflc/(umax*dx1+vmax*dy1+wmax*dz1)
      endif
      dtv=cflv*Re/(dx2+dy2+dz2)

      dt=min(dtc,dtv,dtcfl)
      dt1=1.0d0/dt
      dt2=dt1*dt1

      t=t+dt

      return
      end


c-----------------------------------------------------------------------
