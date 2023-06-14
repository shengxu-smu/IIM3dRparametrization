c-----------------------------------------------------------------------
c
      subroutine qjc_rhs
      include 'parameter.inc'
      include 'surface.inc'
      integer isphr
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,xe1,ye1,ze1,xe2,ye2,ze2
      real*8 gacobi,f1,f2,f3,fn,fo
      real*8 foo(0:ns1,0:ns2),boo(0:ns1,0:ns2)

      DO m=1,ms

      isphr=i01(m)

c  pjc0-fn
      do m2=0,ns2
        do m1=0,ns1
          pjc0(m1,m2,m)=fp(m1,m2,m)
        enddo
      enddo

c  10-fn/a1, 14-fn/a2
      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=fp(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,1,1)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,10,m)=foo(m1,m2)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,14,m)=foo(m1,m2)
        enddo
      enddo

c  15-fn/a1a1
      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=fr(m1,m2,10,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,1,2)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,15,m)=foo(m1,m2)
        enddo
      enddo

c  16-fn/a1a2, 17-fn/a2a2
      do n=2,3
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fr(m1,m2,4*n+2,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,2,1)
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,14+n,m)=foo(m1,m2)
          enddo
        enddo
      enddo

c 20-(\tilde{f1}/a1+\tilde{f2}/a2)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,20,m)=0.0d0
        enddo
      enddo

c  18-(\tilde{f1}/a1+\tilde{f2}/a2)/a1
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,18,m)=0.0d0
          enddo
        enddo

c  19-(\tilde{f1}/a1+\tilde{f2}/a2)/a2
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,19,m)=0.0d0
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
