c-----------------------------------------------------------------------
c
      subroutine time_output
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      real*8 tx,ty,bx,by,sz,f00,f01,f10,f11,fx,fy,fz
      real*8 af(ms),cx(ms),cy(ms),cz(ms),vol(ms)
      real*8 tqx(ms),tqy(ms),tqz(ms)

      DO m=1,ms
        if(i01(m).eq.1) then
          af(m)=pi*ea(m)*eb(m)
        elseif(i01(m).eq.0) then
          af(m)=pi*(ea(m)+ec(m))*(eb(m)+ec(m))-
     .          pi*(ea(m)-ec(m))*(eb(m)-ec(m))
        endif
        cx(m)=0.0d0
        cy(m)=0.0d0
        cz(m)=0.0d0
        vol(m)=0.0d0
        tqx(m)=0.0d0
        tqy(m)=0.0d0
        tqz(m)=0.0d0
        do m2=0,ns2-1
          do m1=0,ns1-1
            f00=-fr(m1,m2,1,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,7,m)
            f01=-fr(m1,m2+1,1,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,7,m)
            f10=-fr(m1+1,m2,1,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,7,m)
            f11=-fr(m1+1,m2+1,1,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,7,m)
            fx=0.25d0*(f00+f01+f10+f11)
            cx(m)=cx(m)-0.25d0*(f00+f01+f10+f11)
            f00=-fr(m1,m2,2,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,8,m)
            f01=-fr(m1,m2+1,2,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,8,m)
            f10=-fr(m1+1,m2,2,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,8,m)
            f11=-fr(m1+1,m2+1,2,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,8,m)
            fy=0.25d0*(f00+f01+f10+f11)
            cy(m)=cy(m)-0.25d0*(f00+f01+f10+f11)
            f00=-fr(m1,m2,3,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,9,m)
            f01=-fr(m1,m2+1,3,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,9,m)
            f10=-fr(m1+1,m2,3,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,9,m)
            f11=-fr(m1+1,m2+1,3,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,9,m)
            fz=0.25d0*(f00+f01+f10+f11)
            cz(m)=cz(m)-0.25d0*(f00+f01+f10+f11)

            tqx(m)=tqx(m)-
     .             ((ys(m1,m2,m)-ysc(m))*fz-(zs(m1,m2,m)-zsc(m))*fy)
            tqy(m)=tqy(m)-
     .             ((zs(m1,m2,m)-zsc(m))*fx-(xs(m1,m2,m)-xsc(m))*fz)
            tqz(m)=tqz(m)-
     .             ((xs(m1,m2,m)-xsc(m))*fy-(ys(m1,m2,m)-ysc(m))*fx)

            tx=xs(m1+1,m2,m)-xs(m1,m2,m)
            ty=ys(m1+1,m2,m)-ys(m1,m2,m)
            bx=xs(m1,m2+1,m)-xs(m1,m2,m)
            by=ys(m1,m2+1,m)-ys(m1,m2,m)
            sz=(zs(m1,m2,m)+zs(m1+1,m2,m)+zs(m1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
            tx=xs(m1+1,m2+1,m)-xs(m1,m2+1,m)
            ty=ys(m1+1,m2+1,m)-ys(m1,m2+1,m)
            bx=xs(m1+1,m2+1,m)-xs(m1+1,m2,m)
            by=ys(m1+1,m2+1,m)-ys(m1+1,m2,m)
            sz=(zs(m1+1,m2,m)+zs(m1+1,m2+1,m)+zs(m1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
           enddo
        enddo
        if(i01(m).eq.1) then
          cx(m)=2.0d0*cx(m)*dalfa11*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa11*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa11*dalfa20/af(m)
          tqx(m)=2.0d0*tqx(m)*dalfa11*dalfa20
          tqy(m)=2.0d0*tqy(m)*dalfa11*dalfa20
          tqz(m)=2.0d0*tqz(m)*dalfa11*dalfa20
          do m2=1,ns2-2
            tx=xs(0,m2+1,m)-xs(0,0,m)
            ty=ys(0,m2+1,m)-ys(0,0,m)
            bx=xs(0,m2+1,m)-xs(0,m2,m)
            by=ys(0,m2+1,m)-ys(0,m2,m)
            sz=(zs(0,0,m)+zs(0,m2,m)+zs(0,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)          
            tx=xs(ns1,0,m)-xs(ns1,m2+1,m)
            ty=ys(ns1,0,m)-ys(ns1,m2+1,m)
            bx=xs(ns1,m2+1,m)-xs(ns1,m2,m)
            by=ys(ns1,m2+1,m)-ys(ns1,m2,m)
            sz=(zs(ns1,0,m)+zs(ns1,m2,m)+zs(ns1,m2+1,m))/3.0d0
            vol(m)=vol(m)+0.5d0*sz*(tx*by-ty*bx)
          enddo
        else
          cx(m)=2.0d0*cx(m)*dalfa10*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa10*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa10*dalfa20/af(m)
          tqx(m)=2.0d0*tqx(m)*dalfa10*dalfa20
          tqy(m)=2.0d0*tqy(m)*dalfa10*dalfa20
          tqz(m)=2.0d0*tqz(m)*dalfa10*dalfa20
        endif
        cx(m)=cx(m)+2.0d0*vol(m)*xsctt(m)/af(m)
        cy(m)=cy(m)+2.0d0*vol(m)*ysctt(m)/af(m)
        cz(m)=cz(m)+2.0d0*vol(m)*zsctt(m)/af(m)
      ENDDO

      write(61,100)t,(cx(m),cy(m),cz(m),tqx(m),tqy(m),tqz(m),m=1,ms)
 100  format(1x,200e16.6e4)

      call fn_correct

      DO m=1,ms
        cx(m)=0.0d0
        cy(m)=0.0d0
        cz(m)=0.0d0
        tqx(m)=0.0d0
        tqy(m)=0.0d0
        tqz(m)=0.0d0
        do m2=0,ns2-1
          do m1=0,ns1-1
            f00=-fr(m1,m2,1,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,7,m)
            f01=-fr(m1,m2+1,1,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,7,m)
            f10=-fr(m1+1,m2,1,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,7,m)
            f11=-fr(m1+1,m2+1,1,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,7,m)
            fx=0.25d0*(f00+f01+f10+f11)
            cx(m)=cx(m)-0.25d0*(f00+f01+f10+f11)
            f00=-fr(m1,m2,2,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,8,m)
            f01=-fr(m1,m2+1,2,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,8,m)
            f10=-fr(m1+1,m2,2,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,8,m)
            f11=-fr(m1+1,m2+1,2,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,8,m)
            fy=0.25d0*(f00+f01+f10+f11)
            cy(m)=cy(m)-0.25d0*(f00+f01+f10+f11)
            f00=-fr(m1,m2,3,m)/Re+fr(m1,m2,4,m)*fst(m1,m2,9,m)
            f01=-fr(m1,m2+1,3,m)/Re+fr(m1,m2+1,4,m)*fst(m1,m2+1,9,m)
            f10=-fr(m1+1,m2,3,m)/Re+fr(m1+1,m2,4,m)*fst(m1+1,m2,9,m)
            f11=-fr(m1+1,m2+1,3,m)/Re+
     .           fr(m1+1,m2+1,4,m)*fst(m1+1,m2+1,9,m)
            fz=0.25d0*(f00+f01+f10+f11)
            cz(m)=cz(m)-0.25d0*(f00+f01+f10+f11)

            tqx(m)=tqx(m)-
     .             ((ys(m1,m2,m)-ysc(m))*fz-(zs(m1,m2,m)-zsc(m))*fy)
            tqy(m)=tqy(m)-
     .             ((zs(m1,m2,m)-zsc(m))*fx-(xs(m1,m2,m)-xsc(m))*fz)
            tqz(m)=tqz(m)-
     .             ((xs(m1,m2,m)-xsc(m))*fy-(ys(m1,m2,m)-ysc(m))*fx)
           enddo
        enddo

        if(i01(m).eq.1) then
          cx(m)=2.0d0*cx(m)*dalfa11*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa11*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa11*dalfa20/af(m)
          tqx(m)=2.0d0*tqx(m)*dalfa11*dalfa20
          tqy(m)=2.0d0*tqy(m)*dalfa11*dalfa20
          tqz(m)=2.0d0*tqz(m)*dalfa11*dalfa20
        else
          cx(m)=2.0d0*cx(m)*dalfa10*dalfa20/af(m)
          cy(m)=2.0d0*cy(m)*dalfa10*dalfa20/af(m)
          cz(m)=2.0d0*cz(m)*dalfa10*dalfa20/af(m)
          tqx(m)=2.0d0*tqx(m)*dalfa10*dalfa20
          tqy(m)=2.0d0*tqy(m)*dalfa10*dalfa20
          tqz(m)=2.0d0*tqz(m)*dalfa10*dalfa20
        endif
      ENDDO

      write(62,100)t,(cx(m),cy(m),cz(m),tqx(m),tqy(m),tqz(m),m=1,ms)

      return
      end


c-----------------------------------------------------------------------
