c-----------------------------------------------------------------------
c
      subroutine interpolate(xa,ya,np,x,y,dy)
      implicit none
      integer np,nmax,ns,n,m,i
      parameter(nmax=10)
      real*8 dy,x,y,xa(np),ya(np)
      real*8 den,dif,dift,ho,hp,w,c(nmax),d(nmax)

      n=np
      if(n.gt.nmax) then
        write(*,*)'  !! erro: interpolation order too high !!'
        stop
      endif

      ns=1
      dif=abs(x-xa(1))

      do i=1,n
        dift=abs(x-xa(i))
        if(dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
      enddo

      y=ya(ns)
      ns=ns-1

      do m=1,n-1
        do i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.0d0) then
            write(*,*)'  !! error: failure in interpolation !!'
            stop
          endif
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
        enddo
        if(2*ns.lt.n-m) then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
      enddo

      return
      end

c-----------------------------------------------------------------------
