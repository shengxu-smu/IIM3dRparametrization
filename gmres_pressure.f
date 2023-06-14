      subroutine gmres_pressure
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
c      include "mkl_rci.fi"
      integer ize,nknown,itercount,irc_request,ks,isphr
      parameter (ize=128,nknown=kks*ms)
      integer ipar(ize)
      real*8 dvar,goo
      real*8 dpar(ize)
      real*8 tmp(nknown*(2*nknown+1)+(nknown*(nknown+9))/2+1)
      real*8 rhs(nknown),rhsb(nknown),solution(nknown),residual(nknown)
      real*8 boo(0:ks1-1,0:ks2-1),foo(0:ns1,0:ns2)
c
c---------------------------------------------------------------------------
c External BLAS functions are taken from MKL BLAS to use
c with the RCI (P)FGMRES solver
c---------------------------------------------------------------------------
      real*8 dnrm2
      external dnrm2

c---------------------------------------------------------------------------
c Initialize the right hand side
c---------------------------------------------------------------------------
      DO m=1,ms
      ks=(m-1)*kks

      do k=1,kks
        m1=mod(k-1,ks1)*kip1
        m2=int((k-1)/ks1)*kip2
        rhs(ks+k)=fq(m1,m2,m)
        solution(ks+k)=fp(m1,m2,m)
      enddo

      ENDDO

c---------------------------------------------------------------------------
c Save the right-hand side in vector rhsb for future use
c---------------------------------------------------------------------------
      call dcopy(nknown,rhs,1,rhsb,1)

c---------------------------------------------------------------------------
c Initialize the solver
c---------------------------------------------------------------------------
      call dfgmres_init(nknown,solution,rhs,irc_request,ipar,dpar,tmp)
      if(irc_request.ne.0) goto 999

c---------------------------------------------------------------------------
c Set the desired parameters:
c do the restart after 5 iterations
c LOGICAL parameters:
c do not do the stopping test for the maximal number of iterations
c do not do the Preconditioned iterations of FGMRES method
c DOUBLE PRECISION parameters
c set the relative tolerance to 1.0d-3 instead of default value 1.0d-6
c---------------------------------------------------------------------------
      ipar(15)=196
      ipar(5)=196
      ipar(8)=1
      ipar(11)=0
      dpar(1)=1.0d-3
c
c---------------------------------------------------------------------------
c Check the correctness and consistency of the newly set parameters
c---------------------------------------------------------------------------
      call dfgmres_check(nknown,solution,rhs,irc_request,ipar,dpar,tmp)
      if(irc_request.ne.0) goto 999

c---------------------------------------------------------------------------
c Compute the solution by RCI (P)FGMRES solver with
c Reverse Communication starts here
c---------------------------------------------------------------------------
1     call dfgmres(nknown,solution,rhs,irc_request,ipar,dpar,tmp)

c---------------------------------------------------------------------------
c If irc_request=0, then the solution was found with the required precision
c---------------------------------------------------------------------------
      if(irc_request.eq.0) goto 3

c---------------------------------------------------------------------------
c If irc_request=1, then compute the matrix-vector prodct A*tmp(ipar(22))
c and put the result in vector tmp(ipar(23))
c---------------------------------------------------------------------------

      if(irc_request.eq.1) then

      DO m=1,ms

      isphr=i01(m)
      ks=(m-1)*kks

      do k2=0,ks2-1
        do k1=0,ks1-1
          boo(k1,k2)=tmp(ipar(22)+ks+k2*ks1+k1)
        enddo
      enddo
      call spline_interpolate(boo,foo,isphr)
      call filter(foo,isphr,32,32,1.0d0)
      goo=0.0d0
      if(isphr.eq.1) then
        do m2=0,ns2-1
          do m1=0,ns1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble((ns1+1)*ns2)
      else
        do m2=0,ns2-1
          do m1=0,ns1-1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble(ns1*ns2)
      endif
      do m2=0,ns2
        do m1=0,ns1
          fp(m1,m2,m)=foo(m1,m2)-goo
        enddo
      enddo

      ENDDO

 200  format(1x,200e16.6)

      call gmres_mvproduct

      DO m=1,ms

      ks=(m-1)*kks

      do k=1,kks
        m1=mod(k-1,ks1)*kip1
        m2=int((k-1)/ks1)*kip2
        tmp(ipar(23)+ks+k-1)=fq(m1,m2,m)
      enddo

      ENDDO

      goto 1

      endif

c---------------------------------------------------------------------------
c If irc_request=2, then do the user-defined stopping test
c The residual stopping test for the computed solution is performed here
c---------------------------------------------------------------------------
c NOTE: from this point vector rhsb is no longer containing the right-hand
c side of the problem! It contains the current FGMRES approximation to the
c solution. If you need to keep the right-hand side, save it in some other
c vector before the call to DFGMRES routine. Here we saved it in vector
c rhs. The vector rhsb is used instead of rhs to preserve the original
c right-hand side of the problem and guarantee the proper restart of FGMRES
c method. Vector rhsb will be altered when computing the residual stopping
c criterion!
c---------------------------------------------------------------------------
      if(irc_request.eq.2) then
c Request to the DFGMRES_GET routine to put the solution into rhsb via ipar(13)
        ipar(13)=1
c Get the current FGMRES solution in the vector rhsb
        call dfgmres_get(nknown,solution,rhsb,irc_request,ipar,
     .                   dpar,tmp,itercount)
c Compute the current true residual
      DO m=1,ms

      isphr=i01(m)
      ks=(m-1)*kks

      do k2=0,ks2-1
        do k1=0,ks1-1
          boo(k1,k2)=rhsb(1+ks+k2*ks1+k1)
        enddo
      enddo
      call spline_interpolate(boo,foo,isphr)
      call filter(foo,isphr,32,32,1.0d0)
      goo=0.0d0
      if(isphr.eq.1) then
        do m2=0,ns2-1
          do m1=0,ns1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble((ns1+1)*ns2)
      else
        do m2=0,ns2-1
          do m1=0,ns1-1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble(ns1*ns2)
      endif
      do m2=0,ns2
        do m1=0,ns1
          fp(m1,m2,m)=foo(m1,m2)-goo
        enddo
      enddo

      ENDDO

      call gmres_mvproduct

      DO m=1,ms
      ks=(m-1)*kks

      do k=1,kks
        m1=mod(k-1,ks1)*kip1
        m2=int((k-1)/ks1)*kip2
        residual(ks+k)=fq(m1,m2,m)
      enddo

      ENDDO

      call daxpy(nknown,-1.0d0,rhs,1,residual,1)
      dvar=dnrm2(nknown,residual,1)
      if(dvar.lt.1.0d-2) then
        goto 3
      else
        goto 1
      endif

      endif

C---------------------------------------------------------------------------
C If irc_request=4, then check if the norm of the next generated vector is
C not zero up to rounding and computational errors. The norm is contained
C in dpar(7) parameter
C---------------------------------------------------------------------------
      if(irc_request.eq.4) then
        if(dpar(7).lt.1.0d-12) then
          goto 3
        else
          goto 1
        endif
      else
c---------------------------------------------------------------------------
c If irc_request=anything else, then DFGMRES subroutine failed
c to compute the solution vector: solution(nknown)
c---------------------------------------------------------------------------
        goto 999
      endif
c
c---------------------------------------------------------------------------
c Reverse Communication ends here
c Get the current iteration number and the FGMRES solution. (DO NOT FORGET to
c call DFGMRES_GET routine as computed_solution is still containing
c the initial guess!). Request to DFGMRES_GET to put the solution into
c vector solution(nknown) via ipar(13)
c---------------------------------------------------------------------------
3     ipar(13)=0
      call dfgmres_get(nknown,solution,rhs,irc_request,ipar,
     .                 dpar,tmp,itercount)
      print *,' Number of GMRES iterations: ',itercount
      goto 1000

999   print *,'The solver has returned the ERROR code ', irc_request

1000  continue

      DO m=1,ms

      isphr=i01(m)
      ks=(m-1)*kks

      do k2=0,ks2-1
        do k1=0,ks1-1
          boo(k1,k2)=solution(1+ks+k2*ks1+k1)
        enddo
      enddo
      call spline_interpolate(boo,foo,isphr)
      call filter(foo,isphr,32,32,1.0d0)
      goo=0.0d0
      if(isphr.eq.1) then
        do m2=0,ns2-1
          do m1=0,ns1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble((ns1+1)*ns2)
      else
        do m2=0,ns2-1
          do m1=0,ns1-1
            goo=goo+foo(m1,m2)
          enddo
        enddo
        goo=goo/dble(ns1*ns2)
      endif
      do m2=0,ns2
        do m1=0,ns1
          fp(m1,m2,m)=foo(m1,m2)-goo
        enddo
      enddo

      ENDDO

      END
