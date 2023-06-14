c----------------------------------------------------------------------
c
      program main
      include 'parameter.inc'
      include 'field.inc'
      include 'surface.inc'
      character*8 day1,day2
      integer*4 time1(3),time2(3)
      real*8 tmp,tstart,tend
      real*4 etime,total,cost(2)

      if(mod(ns2,2).ne.0) then
        write(*,*)'  !! erro: ns2 is not even !!'
        stop
      endif
      if(ms.gt.1) then
        write(*,*)'  !! note: multiple objects !!'
      endif

      call run_input
      call mesh
      call initial

      nstart=0
      if(iread.eq.1) then
        call data_read
        call surface_derivative
      endif
      nstart=nstart+1
      nend=nstart+nstep-1
      tstart=t
      tmp=t

      call date_and_time(day1)
      call itime(time1)

      do n=nstart,nend
        call cfl
        write(*,*)
        write(*,*)'! n = ',n,' t = ',t
        call rk
        if(itout.eq.1.and.mod(n,1).eq.0) then
          call time_output
        endif
        if(ianimation.eq.1) then
          if(tmp.ge.t-dt.and.tmp.lt.t)then
            call animation
            tmp=tmp+tout
          endif
        endif
      enddo
      tend=t

      call date_and_time(day2)
      call itime(time2)

      if(iwrite.eq.1) then
        call data_write
      endif
      if(iplot.eq.1) then
        call surface_plot
        call field_plot
      endif

      total=etime(cost)
      call run_output(tstart,tend,total,day1,time1,day2,time2)

      write(*,*)
      write(*,*)'start @'
      write(*,*)'date: ', day1
      write(*,2000) time1
      write(*,*)'end @'
      write(*,*)'date: ', day2
      write(*,2000) time2
      write(*,*)'user time   = ',cost(1)/3600.0d0,' hours'
      write(*,*)'system time = ',cost(2)/3600.0d0,' hours'
      write(*,*)'total time  = ',total/3600.0d0,' hours'
 2000 format (' time: ',i2.2, ':', i2.2, ':', i2.2 )

      stop
      end


c-----------------------------------------------------------------------
