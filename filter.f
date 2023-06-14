c-----------------------------------------------------------------------
c
      subroutine filter(f,isphr,n1c,n2c,flip)
      include 'parameter.inc'
      include 'surface.inc'
      include 'fft.inc'
      integer isphr,nh,nd,n1c,n2c,nonlinear,nc
      parameter(nd=2*ns1+2,nh=ns2/2) 
      real*8 flip
      real*8 f(0:ns1,0:ns2)
      
      nonlinear=0

      IF(isphr.EQ.0) THEN

      do m2=0,ns2-1
        do m1=0,ns1-1
          r2(m1+1,m2+1)=f(m1,m2)
        enddo
      enddo

      call vrfftf(ns1,ns2,r2,w2,ns1,wsave2)
      do m2=n2c,ns2
        do m1=1,ns1
          r2(m1,m2)=0.0d0
        enddo
      enddo
      call vrfftb(ns1,ns2,r2,w2,ns1,wsave2)
      
      do m2=1,ns2
        do m1=1,ns1
          r1(m2,m1)=r2(m1,m2)
        enddo
      enddo

      call vrfftf(ns2,ns1,r1,w1,ns2,wsave1)
      do m1=n1c,ns1
        do m2=1,ns2
          r1(m2,m1)=0.0d0
        enddo
      enddo
      call vrfftb(ns2,ns1,r1,w1,ns2,wsave1)

      do m2=0,ns2-1
        do m1=0,ns1-1
          f(m1,m2)=r1(m2+1,m1+1)
        enddo
      enddo
      do m1=0,ns1-1
        f(m1,ns2)=f(m1,0)
      enddo
      do m2=0,ns2
        f(ns1,m2)=f(0,m2)
      enddo

      ENDIF

c

      IF(isphr.EQ.1) THEN

      do m2=1,nh
        do m1=1,ns1+1
          r3(m2,m1)=f(m1-1,m2-1)
        enddo
        do m1=ns1+2,nd
          r3(m2,m1)=flip*f(nd-m1,m2+nh-1)
        enddo
      enddo

      call vrfftf(nh,nd,r3,w3,nh,wsave3)
      do m1=n1c,nd
        do m2=1,nh
          r3(m2,m1)=0.0d0
        enddo
      enddo
      call vrfftb(nh,nd,r3,w3,nh,wsave3)

      do m1=1,ns1+1
        do m2=1,nh
          r4(m1,m2)=r3(m2,m1)
        enddo
      enddo
      do m1=1,ns1+1
        do m2=nh+1,ns2
          r4(m1,m2)=flip*r3(m2-nh,nd-m1+1)
        enddo
      enddo

      call vrfftf(ns1+1,ns2,r4,w4,ns1+1,wsave2)
      if(nonlinear.eq.0) then
        do m2=n2c,ns2
          do m1=1,ns1+1
            r4(m1,m2)=0.0d0
          enddo
        enddo
      else
        do m1=1,ns1+1
          nc=4+int(dble(n2c-4)*
     .       (1.0d0-(2.0d0*dble(m1-1)/dble(ns1)-1.0d0)**8.0d0))
          do m2=nc,ns2
             r4(m1,m2)=0.0d0
          enddo
        enddo
      endif
      call vrfftb(ns1+1,ns2,r4,w4,ns1+1,wsave2)

      do m1=0,ns1
        do m2=0,ns2-1
          f(m1,m2)=r4(m1+1,m2+1)
        enddo
      enddo

      do m1=0,ns1
        f(m1,ns2)=f(m1,0)
      enddo

      ENDIF

      return
      end


c-----------------------------------------------------------------------

     
