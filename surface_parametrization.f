c-----------------------------------------------------------------------
c
      subroutine surface_parametrization
      include 'parameter.inc'
      include 'surface.inc'
      integer isphr,jectob,n2,nc
      real*8 angle1,angle2,rr,rx,ry,rz
      real*8 c,s,e,rxx,rxy,rxz,ryx,ryy,ryz,rzx,rzy,rzz
      real*8 aprime,bprime,ds2,angle,adprime,bdprime,ds1,r
      real*8 xr(0:ns2),yr(0:ns2)
      real*8 foo(0:ns1,0:ns2)


c  1: rounded box; 2: ellipsoidal
      jectob=1

      n2=ns2/2

      DO m=1,ms

      isphr=i01(m)

      if(i01(m).eq.1) then
        IF(jectob.eq.1) THEN
        aprime=ea(m)-eb(m)
        bprime=eb(m)
        ds2=(2.0d0*aprime+pi*bprime)/dble(ns2)
        nc=int(0.5d0*aprime/ds2)
        do m2=0,nc
          xr(m2)=dble(m2)*ds2
          yr(m2)=0.5d0*bprime
        enddo
        do m2=nc+1,n2-nc-1
          angle=(2.0d0*dble(m2)*ds2-aprime)/bprime
          xr(m2)=0.5d0*aprime+0.5d0*bprime*dsin(angle)
          yr(m2)=0.5d0*bprime*dcos(angle)
        enddo
        do m2=n2-nc,n2+nc
          xr(m2)=dble(n2-m2)*ds2
          yr(m2)=-0.5d0*bprime
        enddo
        do m2=n2+nc+1,ns2-nc-1
          angle=(3.0d0*aprime+2.0d0*(pi*bprime-dble(m2)*ds2))/bprime
          xr(m2)=-0.5d0*aprime-0.5d0*bprime*dsin(angle)
          yr(m2)=0.5d0*bprime*dcos(angle)
        enddo
        do m2=ns2-nc,ns2
          xr(m2)=dble(m2-ns2)*ds2
          yr(m2)=0.5d0*bprime
        enddo
        do m2=0,ns2
          adprime=2.0d0*dsqrt(xr(m2)*xr(m2)+yr(m2)*yr(m2))-ec(m)
          bdprime=ec(m)
          ds1=(adprime+0.5d0*pi*bdprime)/dble(ns1+1)
          nc=int((adprime-ds1)/(2.0d0*ds1))
          do m1=0,nc
            r=0.5d0*ds1+dble(m1)*ds1
            xss(m1,ns2-m2,m)=2.0d0*r*xr(m2)/(adprime+ec(m))
            yss(m1,ns2-m2,m)=2.0d0*r*yr(m2)/(adprime+ec(m))
            zss(m1,ns2-m2,m)=0.5d0*bdprime
          enddo
          do m1=nc+1,ns1-nc-1
            angle=(ds1+2.0d0*dble(m1)*ds1-adprime)/bdprime
            r=0.5d0*adprime+0.5d0*bdprime*dsin(angle)
            xss(m1,ns2-m2,m)=2.0d0*r*xr(m2)/(adprime+ec(m))
            yss(m1,ns2-m2,m)=2.0d0*r*yr(m2)/(adprime+ec(m))
            zss(m1,ns2-m2,m)=0.5d0*bdprime*dcos(angle)
          enddo
          do m1=ns1-nc,ns1
            r=0.5d0*ds1+dble(ns1-m1)*ds1
            xss(m1,ns2-m2,m)=2.0d0*r*xr(m2)/(adprime+ec(m))
            yss(m1,ns2-m2,m)=2.0d0*r*yr(m2)/(adprime+ec(m))
            zss(m1,ns2-m2,m)=-0.5d0*bdprime
          enddo
        enddo
        ENDIF
        IF(jectob.eq.2) THEN
        do m2=0,ns2
          angle2=alfa20(m2)
          do m1=0,ns1        
            angle1=alfa11(m1)
            xss(m1,m2,m)=ea(m)*dsin(angle1)*dcos(angle2)
            yss(m1,m2,m)=eb(m)*dsin(angle1)*dsin(angle2)
            zss(m1,m2,m)=ec(m)*dcos(angle1)
          enddo
        enddo
        ENDIF
      else
        do m2=0,ns2
          angle2=alfa20(m2)
          do m1=0,ns1
            angle1=alfa10(m1)
            xss(m1,m2,m)=ec(m)*dcos(angle1)
            yss(m1,m2,m)=(ec(m)*dsin(angle1)+ea(m))*dcos(angle2)
            zss(m1,m2,m)=(ec(m)*dsin(angle1)+eb(m))*dsin(angle2)
          enddo
        enddo
      endif

      do m1=0,ns1
        xss(m1,ns2,m)=xss(m1,0,m)
        yss(m1,ns2,m)=yss(m1,0,m)
        zss(m1,ns2,m)=zss(m1,0,m)
      enddo
      if(i01(m).eq.0) then
        do m2=0,ns2
          xss(ns1,m2,m)=xss(0,m2,m)
          yss(ns1,m2,m)=yss(0,m2,m)
          zss(ns1,m2,m)=zss(0,m2,m)
        enddo
      endif

      do m2=0,ns2
        do m1=0,ns1
          foo(m1,m2)=xss(m1,m2,m)
        enddo
      enddo
      call filter(foo,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          xss(m1,m2,m)=foo(m1,m2)
        enddo
      enddo
      do m2=0,ns2
        do m1=0,ns1
          foo(m1,m2)=yss(m1,m2,m)
        enddo
      enddo
      call filter(foo,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          yss(m1,m2,m)=foo(m1,m2)
        enddo
      enddo
      do m2=0,ns2
        do m1=0,ns1
          foo(m1,m2)=zss(m1,m2,m)
        enddo
      enddo
      call filter(foo,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          zss(m1,m2,m)=foo(m1,m2)
        enddo
      enddo

      ENDDO

100   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
