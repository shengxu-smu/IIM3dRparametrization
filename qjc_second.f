c-----------------------------------------------------------------------
c
      subroutine qjc_second
      include 'parameter.inc'
      include 'surface.inc'
      integer id(7),ix(6)
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,gacobi,gacobi2,g1,g2
      real*8 xe1,ye1,ze1,xe2,ye2,ze2
      real*8 fo1,fo2,fo3,bo1,bo2,bo3
      real*8 bb(7),dd(7),xx(6)

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
          g1=dsqrt(xt*xt+yt*yt+zt*zt)
          xt=xt/g1
          yt=yt/g1
          zt=zt/g1
          g2=dsqrt(xb*xb+yb*yb+zb*zb)
          xb=xb/g2
          yb=yb/g2
          zb=zb/g2
          bb(1)=fr(m1,m2,15,m)-
     .             (snd(m1,m2,1,m)*pjc1(m1,m2,1,m)+
     .              snd(m1,m2,2,m)*pjc1(m1,m2,2,m)+
     .              snd(m1,m2,3,m)*pjc1(m1,m2,3,m))
          bb(2)=fr(m1,m2,16,m)-
     .             (snd(m1,m2,7,m)*pjc1(m1,m2,1,m)+
     .              snd(m1,m2,8,m)*pjc1(m1,m2,2,m)+
     .              snd(m1,m2,9,m)*pjc1(m1,m2,3,m))
          bb(3)=fr(m1,m2,18,m)-
     .             (snd(m1,m2,10,m)*pjc1(m1,m2,1,m)+
     .              snd(m1,m2,11,m)*pjc1(m1,m2,2,m)+
     .              snd(m1,m2,12,m)*pjc1(m1,m2,3,m))
          bb(4)=bb(2)
          bb(5)=fr(m1,m2,17,m)-
     .             (snd(m1,m2,4,m)*pjc1(m1,m2,1,m)+
     .              snd(m1,m2,5,m)*pjc1(m1,m2,2,m)+
     .              snd(m1,m2,6,m)*pjc1(m1,m2,3,m))
          bb(6)=fr(m1,m2,19,m)-
     .             (snd(m1,m2,13,m)*pjc1(m1,m2,1,m)+
     .              snd(m1,m2,14,m)*pjc1(m1,m2,2,m)+
     .              snd(m1,m2,15,m)*pjc1(m1,m2,3,m))
          bb(7)=0.0d0
          fo1=(xe1*bb(1)+xe2*bb(2)+xn*bb(3))*gacobi2
          fo2=(ye1*bb(1)+ye2*bb(2)+yn*bb(3))*gacobi2
          fo3=(ze1*bb(1)+ze2*bb(2)+zn*bb(3))*gacobi2
          bb(1)=fo1/g1
          bb(2)=fo2/g1
          bb(3)=fo3/g1
          fo1=(xe1*bb(4)+xe2*bb(5)+xn*bb(6))*gacobi2
          fo2=(ye1*bb(4)+ye2*bb(5)+yn*bb(6))*gacobi2
          fo3=(ze1*bb(4)+ze2*bb(5)+zn*bb(6))*gacobi2
          bb(4)=fo1/g2
          bb(5)=fo2/g2
          bb(6)=fo3/g2
          if(abs(xt).ge.abs(yt).and.abs(xt).ge.abs(zt)) then
            id(1)=7
            id(2)=1
            id(3)=4
            id(4)=2
            id(5)=5
            id(6)=3
            id(7)=6
            ix(1)=1
            ix(2)=2
            ix(3)=3
            ix(4)=4
            ix(5)=5
            ix(6)=6
            fo1=xt
            fo2=yt
            fo3=zt
            bo1=xb
            bo2=yb
            bo3=zb
          endif
          if(abs(yt).ge.abs(xt).and.abs(yt).ge.abs(zt)) then
            id(1)=7
            id(2)=2
            id(3)=5
            id(4)=3
            id(5)=6
            id(6)=1
            id(7)=4
            ix(1)=6
            ix(2)=3
            ix(3)=5
            ix(4)=1
            ix(5)=2
            ix(6)=4
            fo1=yt
            fo2=zt
            fo3=xt
            bo1=yb
            bo2=zb
            bo3=xb
          endif
          if(abs(zt).ge.abs(xt).and.abs(zt).ge.abs(yt)) then
            id(1)=7
            id(2)=3
            id(3)=6
            id(4)=2
            id(5)=5
            id(6)=1
            id(7)=4
            ix(1)=6
            ix(2)=5
            ix(3)=3
            ix(4)=4
            ix(5)=2
            ix(6)=1
            fo1=zt
            fo2=yt
            fo3=xt
            bo1=zb
            bo2=yb
            bo3=xb
          endif
          do i=1,7
            dd(i)=bb(id(i))
          enddo
          call c2_solver(fo1,fo2,fo3,bo1,bo2,bo3,dd,xx,0)
          do i=1,6
            pjc2(m1,m2,i,m)=xx(ix(i))
          enddo
          do i=1,6
            pjc(m1,m2,i,m)=0.0d0
          enddo
        enddo
      enddo

      ENDDO

      return
      end

c-----------------------------------------------------------------------
