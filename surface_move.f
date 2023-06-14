c-----------------------------------------------------------------------
c
      subroutine surface_move
      include 'parameter.inc'
      include 'surface.inc'
      real*8 si,ci,sii,cii,siii,ciii,s1,c1,s2,c2,s3,c3
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,gacobi
      real*8 xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3
      real*8 xx(5),yy(5),zz(5)
      real*8 x1(5),y1(5),z1(5),x2(5),y2(5),z2(5),x3(5),y3(5),z3(5)

      t0=t0+0.5d0*dt

      DO m=1,ms

      si=dsin(phi(m))
      ci=dcos(phi(m))
      sii=dsin(theta(m))
      cii=dcos(theta(m))
      siii=dsin(psi(m))
      ciii=dcos(psi(m))

      xsc(m)=xsc0(m)
      ysc(m)=ysc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dcos(pi/3.0d0)
      zsc(m)=zsc0(m)+1.25d0*(dcos(0.8d0*t0)+1.0d0)*dsin(pi/3.0d0)
      xsct(m)=0.0d0
      ysct(m)=0.0d0+1.25d0*(-0.8d0*dsin(0.8d0*t0))*dcos(pi/3.0d0)
      zsct(m)=0.0d0+1.25d0*(-0.8d0*dsin(0.8d0*t0))*dsin(pi/3.0d0)
      xsctt(m)=0.0d0
      ysctt(m)=0.0d0+1.25d0*(-0.64d0*dcos(0.8d0*t0))*dcos(pi/3.0d0)
      zsctt(m)=0.0d0+1.25d0*(-0.64d0*dcos(0.8d0*t0))*dsin(pi/3.0d0)

      phi(m)=phi0(m)
      theta(m)=theta0(m)+pi-0.25d0*pi*(1.0d0-dsin(0.8d0*t0))
      psi(m)=psi0(m)
      phit(m)=0.0d0
      thetat(m)=0.0d0+0.25d0*pi*(0.8d0*dcos(0.8d0*t0))
      psit(m)=0.0d0
      phitt(m)=0.0d0
      thetatt(m)=0.0d0+0.25d0*pi*(-0.64d0*dsin(0.8d0*t0))
      psitt(m)=0.0d0

      s1=dsin(phi(m))
      c1=dcos(phi(m))
      s2=dsin(theta(m))
      c2=dcos(theta(m))
      s3=dsin(psi(m))
      c3=dcos(psi(m))

      omegax(m)=s3*s2*phit(m)+c3*thetat(m)
      omegay(m)=s3*thetat(m)-c3*s2*phit(m)
      omegaz(m)=c2*phit(m)+psit(m)
      omegaxt(m)=s3*s2*phitt(m)+s3*c2*thetat(m)*phit(m)+
     .           c3*psit(m)*s2*phit(m)+
     .           c3*thetatt(m)-s3*psit(m)*thetat(m)
      omegayt(m)=s3*thetatt(m)+c3*psit(m)*thetat(m)
     .           -c3*s2*phitt(m)-c3*c2*thetat(m)*phit(m)
     .           +s3*psit(m)*phit(m)
      omegazt(m)=c2*phitt(m)-s2*thetat(m)*phit(m)+psitt(m)
      
      do m2=0,ns2
        do m1=0,ns1
          xx1=c1*xss(m1,m2,m)-s1*yss(m1,m2,m)
          yy1=s1*xss(m1,m2,m)+c1*yss(m1,m2,m)
          zz1=zss(m1,m2,m)
          xx2=xx1
          yy2=c2*yy1-s2*zz1
          zz2=s2*yy1+c2*zz1
          xx3=c3*xx2-s3*yy2
          yy3=s3*xx2+c3*yy2
          zz3=zz2
          xs(m1,m2,m)=xsc(m)+xx3
          ys(m1,m2,m)=ysc(m)+yy3
          zs(m1,m2,m)=zsc(m)+zz3
          us(m1,m2,m)=xsct(m)+omegay(m)*zz3-omegaz(m)*yy3
          vs(m1,m2,m)=ysct(m)+omegaz(m)*xx3-omegax(m)*zz3
          ws(m1,m2,m)=zsct(m)+omegax(m)*yy3-omegay(m)*xx3

          x3(1)=fst(m1,m2,1,m)
          y3(1)=fst(m1,m2,2,m)
          z3(1)=fst(m1,m2,3,m)
          x3(2)=fst(m1,m2,4,m)
          y3(2)=fst(m1,m2,5,m)
          z3(2)=fst(m1,m2,6,m)
          x3(3)=snd(m1,m2,1,m)
          y3(3)=snd(m1,m2,2,m)
          z3(3)=snd(m1,m2,3,m)
          x3(4)=snd(m1,m2,4,m)
          y3(4)=snd(m1,m2,5,m)
          z3(4)=snd(m1,m2,6,m)
          x3(5)=snd(m1,m2,7,m)
          y3(5)=snd(m1,m2,8,m)
          z3(5)=snd(m1,m2,9,m)

          do n=1,5
            x2(n)=ciii*x3(n)+siii*y3(n)
            y2(n)=-siii*x3(n)+ciii*y3(n)
            z2(n)=z3(n)
            x1(n)=x2(n)
            y1(n)=cii*y2(n)+sii*z2(n)
            z1(n)=-sii*y2(n)+cii*z2(n)
            xx(n)=ci*x1(n)+si*y1(n)
            yy(n)=-si*x1(n)+ci*y1(n)
            zz(n)=z1(n)

            x1(n)=c1*xx(n)-s1*yy(n)
            y1(n)=s1*xx(n)+c1*yy(n)
            z1(n)=zz(n)
            x2(n)=x1(n)
            y2(n)=c2*y1(n)-s2*z1(n)
            z2(n)=s2*y1(n)+c2*z1(n)
            x3(n)=c3*x2(n)-s3*y2(n)
            y3(n)=s3*x2(n)+c3*y2(n)
            z3(n)=z2(n)
          enddo

          fst(m1,m2,1,m)=x3(1)
          fst(m1,m2,2,m)=y3(1)
          fst(m1,m2,3,m)=z3(1)
          fst(m1,m2,4,m)=x3(2)
          fst(m1,m2,5,m)=y3(2)
          fst(m1,m2,6,m)=z3(2)
          snd(m1,m2,1,m)=x3(3)
          snd(m1,m2,2,m)=y3(3)
          snd(m1,m2,3,m)=z3(3)
          snd(m1,m2,4,m)=x3(4)
          snd(m1,m2,5,m)=y3(4)
          snd(m1,m2,6,m)=z3(4)
          snd(m1,m2,7,m)=x3(5)
          snd(m1,m2,8,m)=y3(5)
          snd(m1,m2,9,m)=z3(5)

          xt=x3(1)
          yt=y3(1)
          zt=z3(1)
          xb=x3(2)
          yb=y3(2)
          zb=z3(2)
          xn=yt*zb-zt*yb
          yn=zt*xb-xt*zb
          zn=xt*yb-yt*xb
          gacobi=dsqrt(xn*xn+yn*yn+zn*zn)
          if(gacobi.eq.0.0d0) then
            write(*,*)'  !! erro: pole singularity !!'
            stop
          else
            fst(m1,m2,7,m)=xn
            fst(m1,m2,8,m)=yn
            fst(m1,m2,9,m)=zn
            fst(m1,m2,10,m)=gacobi
          endif

          snd(m1,m2,10,m)=yt*z3(5)+zb*y3(3)-
     .                    zt*y3(5)-yb*z3(3)
          snd(m1,m2,11,m)=zt*x3(5)+xb*z3(3)-
     .                    xt*z3(5)-zb*x3(3)
          snd(m1,m2,12,m)=xt*y3(5)+yb*x3(3)-
     .                    yt*x3(5)-xb*y3(3)
          snd(m1,m2,13,m)=yt*z3(4)+zb*y3(5)-
     .                    zt*y3(4)-yb*z3(5)
          snd(m1,m2,14,m)=zt*x3(4)+xb*z3(5)-
     .                    xt*z3(4)-zb*x3(5)
          snd(m1,m2,15,m)=xt*y3(4)+yb*x3(5)-
     .                    yt*x3(4)-xb*y3(5)
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
