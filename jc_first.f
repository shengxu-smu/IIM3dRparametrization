c-----------------------------------------------------------------------
c
      subroutine jc_first
      include 'parameter.inc'
      include 'surface.inc'
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,xe1,ye1,ze1,xe2,ye2,ze2
      real*8 f1,f2,f3,fo,gacobi,gacobi2

      DO m=1,ms

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
          gacobi2=1.0d0/(gacobi*gacobi)
          xe1=(yb*zn-zb*yn)
          ye1=(zb*xn-xb*zn)
          ze1=(xb*yn-yb*xn)
          xe2=(yn*zt-zn*yt)
          ye2=(zn*xt-xn*zt)
          ze2=(xn*yt-yn*xt)
          f3=fr(m1,m2,1,m)
          ujc1(m1,m2,1,m)=xn*f3*gacobi2
          ujc1(m1,m2,2,m)=yn*f3*gacobi2
          ujc1(m1,m2,3,m)=zn*f3*gacobi2
          f3=fr(m1,m2,2,m)
          vjc1(m1,m2,1,m)=xn*f3*gacobi2
          vjc1(m1,m2,2,m)=yn*f3*gacobi2
          vjc1(m1,m2,3,m)=zn*f3*gacobi2
          f3=fr(m1,m2,3,m)
          wjc1(m1,m2,1,m)=xn*f3*gacobi2
          wjc1(m1,m2,2,m)=yn*f3*gacobi2
          wjc1(m1,m2,3,m)=zn*f3*gacobi2
          f1=fr(m1,m2,10,m)
          f2=fr(m1,m2,14,m)
          f3=fr(m1,m2,20,m)
          pjc1(m1,m2,1,m)=(xe1*f1+xe2*f2+xn*f3)*gacobi2
          pjc1(m1,m2,2,m)=(ye1*f1+ye2*f2+yn*f3)*gacobi2
          pjc1(m1,m2,3,m)=(ze1*f1+ze2*f2+zn*f3)*gacobi2
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
