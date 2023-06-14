c-----------------------------------------------------------------------
c
      subroutine jc_rhs
      include 'parameter.inc'
      include 'surface.inc'
      integer isphr
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,xe1,ye1,ze1,xe2,ye2,ze2
      real*8 gacobi,f1,f2,f3,fn,fo
      real*8 foo(0:ns1,0:ns2),boo(0:ns1,0:ns2)

      DO m=1,ms

      isphr=i01(m)

c  1-ft1, 2-ft2, 3-ft3, 4-fn, 5-\tilde{f1}, 6-\tilde{f2}
      do m2=0,ns2
        do m1=0,ns1
          xt=fst(m1,m2,1,m)
          yt=fst(m1,m2,2,m)
          zt=fst(m1,m2,3,m)
          xb=fst(m1,m2,4,m)
          yb=fst(m1,m2,5,m)
          zb=fst(m1,m2,6,m)
          xn=fst(m1,m2,7,m)
          yn=fst(m1,m2,8,m)
          zn=fst(m1,m2,9,m)
          gacobi=fst(m1,m2,10,m)
          xe1=yb*zn-zb*yn
          ye1=zb*xn-xb*zn
          ze1=xb*yn-yb*xn
          xe2=yn*zt-zn*yt
          ye2=zn*xt-xn*zt
          ze2=xn*yt-yn*xt
          pjc0(m1,m2,m)=fr(m1,m2,4,m)
          f1=-fr(m1,m2,1,m)/(Re*gacobi)
          f2=-fr(m1,m2,2,m)/(Re*gacobi)
          f3=-fr(m1,m2,3,m)/(Re*gacobi)
          fr(m1,m2,5,m)=(f1*xe1+f2*ye1+f3*ze1)/gacobi
          fr(m1,m2,6,m)=(f1*xe2+f2*ye2+f3*ze2)/gacobi
        enddo
      enddo

c  7-ft1/a1, 8-ft2/a1, 9-ft3/a1, 10-fn/a1
c  11-ft1/a2, 12-ft2/a2, 13-ft3/a2, 14-fn/a2
      do n=1,4
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fr(m1,m2,n,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,1,1)
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,n+6,m)=foo(m1,m2)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,2,1)
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,n+10,m)=foo(m1,m2)
          enddo
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

c  18-\tilde{f1}/a1
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fr(m1,m2,5,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,1,1)
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,18,m)=foo(m1,m2)
          enddo
        enddo

c 19-\tilde{f2}/a2, 20-(\tilde{f1}/a1+\tilde{f2}/a2)
      do m2=0,ns2
         do m1=0,ns1
          boo(m1,m2)=fr(m1,m2,6,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fo=foo(m1,m2)
          fr(m1,m2,19,m)=fo
          fr(m1,m2,20,m)=fr(m1,m2,18,m)+fo
        enddo
      enddo

c  18-(\tilde{f1}/a1+\tilde{f2}/a2)/a1
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fr(m1,m2,20,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,1,1)
        do m2=0,ns2
          do m1=0,ns1
            fr(m1,m2,18,m)=foo(m1,m2)
          enddo
        enddo

c  19-(\tilde{f1}/a1+\tilde{f2}/a2)/a2
      do m2=0,ns2
         do m1=0,ns1
          boo(m1,m2)=fr(m1,m2,20,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,19,m)=foo(m1,m2)
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
