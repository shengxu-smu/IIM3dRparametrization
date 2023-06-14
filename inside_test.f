c-----------------------------------------------------------------------
c
      subroutine inside_test(aa,bb,cc,dd,xx,yy,zz,inside)
      implicit none
      integer inside,ia,ib,ic,id
      real*8 xx,yy,zz,xt,yt,zt,xb,yb,zb,xn,yn,zn,side
      real*8 aa(3),bb(3),cc(3),dd(3)

      xt=bb(1)-dd(1)
      yt=bb(2)-dd(2)
      zt=bb(3)-dd(3)
      xb=cc(1)-dd(1)
      yb=cc(2)-dd(2)
      zb=cc(3)-dd(3)
      xn=yt*zb-zt*yb
      yn=zt*xb-xt*zb
      zn=xt*yb-yt*xb
      side=(xn*(xx-dd(1))+yn*(yy-dd(2))+zn*(zz-dd(3)))*
     .     (xn*(aa(1)-dd(1))+yn*(aa(2)-dd(2))+zn*(aa(3)-dd(3)))
      ia=0
      if(side.ge.0.0d0) ia=1

      xt=cc(1)-dd(1)
      yt=cc(2)-dd(2)
      zt=cc(3)-dd(3)
      xb=aa(1)-dd(1)
      yb=aa(2)-dd(2)
      zb=aa(3)-dd(3)
      xn=yt*zb-zt*yb
      yn=zt*xb-xt*zb
      zn=xt*yb-yt*xb
      side=(xn*(xx-dd(1))+yn*(yy-dd(2))+zn*(zz-dd(3)))*
     .     (xn*(bb(1)-dd(1))+yn*(bb(2)-dd(2))+zn*(bb(3)-dd(3)))
      ib=0
      if(side.ge.0.0d0) ib=1

      xt=aa(1)-dd(1)
      yt=aa(2)-dd(2)
      zt=aa(3)-dd(3)
      xb=bb(1)-dd(1)
      yb=bb(2)-dd(2)
      zb=bb(3)-dd(3)
      xn=yt*zb-zt*yb
      yn=zt*xb-xt*zb
      zn=xt*yb-yt*xb
      side=(xn*(xx-dd(1))+yn*(yy-dd(2))+zn*(zz-dd(3)))*
     .     (xn*(cc(1)-dd(1))+yn*(cc(2)-dd(2))+zn*(cc(3)-dd(3)))
      ic=0
      if(side.ge.0.0d0) ic=1

      xt=bb(1)-aa(1)
      yt=bb(2)-aa(2)
      zt=bb(3)-aa(3)
      xb=cc(1)-aa(1)
      yb=cc(2)-aa(2)
      zb=cc(3)-aa(3)
      xn=yt*zb-zt*yb
      yn=zt*xb-xt*zb
      zn=xt*yb-yt*xb
      side=(xn*(xx-aa(1))+yn*(yy-aa(2))+zn*(zz-aa(3)))*
     .     (xn*(dd(1)-aa(1))+yn*(dd(2)-aa(2))+zn*(dd(3)-aa(3)))
      id=0
      if(side.ge.0.0d0) id=1

      inside=0
      if(ia.eq.1.and.ib.eq.1.and.ic.eq.1.and.id.eq.1) inside=1
 
      return
      end


c-----------------------------------------------------------------------

