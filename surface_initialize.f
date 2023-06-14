c-----------------------------------------------------------------------
c
      subroutine surface_initialize
      include 'parameter.inc'
      include 'surface.inc'
      real*8 s1,c1,s2,c2,s3,c3
      real*8 xx1,yy1,zz1,xx2,yy2,zz2,xx3,yy3,zz3

      DO m=1,ms

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
        enddo
      enddo

      ENDDO

      call surface_derivative

      return
      end


c-----------------------------------------------------------------------
