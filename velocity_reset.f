c-----------------------------------------------------------------------
c
      subroutine velocity_reset
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer isphr,imax,imin,jmax,jmin,kmax,kmin,ii,jj,kk,mm
      integer l1,l2,l3,l4,la,lb,lc,ld
      integer iu,id,ju,jd,ku,kd,ie,ic,je,jc,ke,kc
      real*8 s1,c1,s2,c2,s3,c3,xx,yy,zz,x1,y1,z1,x2,y2,z2,x3,y3,z3
      real*8 sn,gacobi,xn,yn,zn,side1,side2,side3,side4,side5
      real*8 xn1,xn2,xn3,xn4,xn5,yn1,yn2,yn3,yn4,yn5,zn1,zn2,zn3,zn4,zn5
      real*8 f1(3),f2(3),f3(3),f4(3),aa(3),bb(3),cc(3),dd(3)

      sn=1.01d0*dsqrt(dx*dx+dy*dy+dz*dz)

      DO m=1,ms

      isphr=i01(m)

      IF(isphr.eq.1) THEN
      imax=min0(int((xsc(m)+ea(m)+dx-x0)/dx)+1,nx)
      imin=max0(int((xsc(m)-ea(m)-dx-x0)/dx)+1,1)

      jmax=min0(int((ysc(m)+ea(m)+dy-y0)/dy)+1,ny)
      jmin=max0(int((ysc(m)-ea(m)-dy-y0)/dy)+1,1)

      kmax=min0(int((zsc(m)+ea(m)+dz-z0)/dz)+1,nz)
      kmin=max0(int((zsc(m)-ea(m)-dz-z0)/dz)+1,1)
      ENDIF

      IF(isphr.eq.0) THEN
      imax=min0(int((xsc(m)+ec(m)+dx-x0)/dx)+1,nx)
      imin=max0(int((xsc(m)-ec(m)-dx-x0)/dx)+1,1)

      jmax=min0(int((ysc(m)+ea(m)+ec(m)+dy-y0)/dy)+1,ny)
      jmin=max0(int((ysc(m)-ea(m)-ec(m)-dy-y0)/dy)+1,1)

      kmax=min0(int((zsc(m)+eb(m)+ec(m)+dz-z0)/dz)+1,nz)
      kmin=max0(int((zsc(m)-eb(m)-ec(m)-dz-z0)/dz)+1,1)
      ENDIF


      s1=dsin(phi(m))
      c1=dcos(phi(m))
      s2=dsin(theta(m))
      c2=dcos(theta(m))
      s3=dsin(psi(m))
      c3=dcos(psi(m))

      do k=kmin,kmax
        do j=jmin,jmax
          do i=imin,imax
            x3=xe(i)-xsc(m)
            y3=yc(j)-ysc(m)
            z3=zc(k)-zsc(m)
            x2=c3*x3+s3*y3
            y2=-s3*x3+c3*y3
            z2=z3
            x1=x2
            y1=c2*y2+s2*z2
            z1=-s2*y2+c2*z2
            xx=c1*x1+s1*y1
            yy=-s1*x1+c1*y1
            zz=z1
            ii=int((xx-x0)/hdx)
            jj=int((yy-y0)/hdy)
            kk=int((zz-z0)/hdz)
            if(ii.ge.0.and.ii.le.2*nx-3) then
              if(jj.ge.0.and.jj.le.2*ny-3) then
                if(kk.ge.0.and.kk.le.2*nz-3) then
                  mm=ind(ii,jj,kk)
                  if(mm.gt.0) then
                    u(i,j,k)=xsct(mm)+omegay(mm)*(zc(k)-zsc(mm))
     .                               -omegaz(mm)*(yc(j)-ysc(mm))
                  endif
                endif
              endif
            endif

            x3=xc(i)-xsc(m)
            y3=ye(j)-ysc(m)
            z3=zc(k)-zsc(m)
            x2=c3*x3+s3*y3
            y2=-s3*x3+c3*y3
            z2=z3
            x1=x2
            y1=c2*y2+s2*z2
            z1=-s2*y2+c2*z2
            xx=c1*x1+s1*y1
            yy=-s1*x1+c1*y1
            zz=z1
            ii=int((xx-x0)/hdx)
            jj=int((yy-y0)/hdy)
            kk=int((zz-z0)/hdz)
            if(ii.ge.0.and.ii.le.2*nx-3) then
              if(jj.ge.0.and.jj.le.2*ny-3) then
                if(kk.ge.0.and.kk.le.2*nz-3) then
                  mm=ind(ii,jj,kk)
                  if(mm.gt.0) then
                    v(i,j,k)=ysct(mm)+omegaz(mm)*(xc(i)-xsc(mm))
     .                               -omegax(mm)*(zc(k)-zsc(mm))
                  endif
                endif
              endif
            endif

            x3=xc(i)-xsc(m)
            y3=yc(j)-ysc(m)
            z3=ze(k)-zsc(m)
            x2=c3*x3+s3*y3
            y2=-s3*x3+c3*y3
            z2=z3
            x1=x2
            y1=c2*y2+s2*z2
            z1=-s2*y2+c2*z2
            xx=c1*x1+s1*y1
            yy=-s1*x1+c1*y1
            zz=z1
            ii=int((xx-x0)/hdx)
            jj=int((yy-y0)/hdy)
            kk=int((zz-z0)/hdz)
            if(ii.ge.0.and.ii.le.2*nx-3) then
              if(jj.ge.0.and.jj.le.2*ny-3) then
                if(kk.ge.0.and.kk.le.2*nz-3) then
                  mm=ind(ii,jj,kk)
                  if(mm.gt.0) then
                    w(i,j,k)=zsct(mm)+omegax(mm)*(yc(j)-ysc(mm))
     .                               -omegay(mm)*(xc(i)-xsc(mm))
                  endif
                endif
              endif
            endif

          enddo
        enddo
      enddo

      ENDDO

c

      DO m=1,ms

      do m2=0,ns2-1
        do m1=0,ns1-1

          f1(1)=xs(m1,m2,m)
          f2(1)=xs(m1+1,m2,m)
          f3(1)=xs(m1+1,m2+1,m)
          f4(1)=xs(m1,m2+1,m)

          f1(2)=ys(m1,m2,m)
          f2(2)=ys(m1+1,m2,m)
          f3(2)=ys(m1+1,m2+1,m)
          f4(2)=ys(m1,m2+1,m)

          f1(3)=zs(m1,m2,m)
          f2(3)=zs(m1+1,m2,m)
          f3(3)=zs(m1+1,m2+1,m)
          f4(3)=zs(m1,m2+1,m)

          gacobi=fst(m1,m2,10,m)
          xn=-fst(m1,m2,7,m)/gacobi
          yn=-fst(m1,m2,8,m)/gacobi
          zn=-fst(m1,m2,9,m)/gacobi
          aa(1)=f1(1)+sn*xn
          aa(2)=f1(2)+sn*yn
          aa(3)=f1(3)+sn*zn

          gacobi=fst(m1+1,m2,10,m)
          xn=-fst(m1+1,m2,7,m)/gacobi
          yn=-fst(m1+1,m2,8,m)/gacobi
          zn=-fst(m1+1,m2,9,m)/gacobi
          bb(1)=f2(1)+sn*xn
          bb(2)=f2(2)+sn*yn
          bb(3)=f2(3)+sn*zn

          gacobi=fst(m1+1,m2+1,10,m)
          xn=-fst(m1+1,m2+1,7,m)/gacobi
          yn=-fst(m1+1,m2+1,8,m)/gacobi
          zn=-fst(m1+1,m2+1,9,m)/gacobi
          cc(1)=f3(1)+sn*xn
          cc(2)=f3(2)+sn*yn
          cc(3)=f3(3)+sn*zn

          gacobi=fst(m1,m2+1,10,m)
          xn=-fst(m1,m2+1,7,m)/gacobi
          yn=-fst(m1,m2+1,8,m)/gacobi
          zn=-fst(m1,m2+1,9,m)/gacobi
          dd(1)=f4(1)+sn*xn
          dd(2)=f4(2)+sn*yn
          dd(3)=f4(3)+sn*zn

c  A-1, B-2, D-4:

          la=int((aa(1)-x0)/hdx)
          lb=int((bb(1)-x0)/hdx)
          ld=int((dd(1)-x0)/hdx)
          l1=int((f1(1)-x0)/hdx)
          l2=int((f2(1)-x0)/hdx)
          l4=int((f4(1)-x0)/hdx)
          iu=max0(la,lb,ld,l1,l2,l4)
          id=min0(la,lb,ld,l1,l2,l4)

          la=int((aa(2)-y0)/hdy)
          lb=int((bb(2)-y0)/hdy)
          ld=int((dd(2)-y0)/hdy)
          l1=int((f1(2)-y0)/hdy)
          l2=int((f2(2)-y0)/hdy)
          l4=int((f4(2)-y0)/hdy)
          ju=max0(la,lb,ld,l1,l2,l4)
          jd=min0(la,lb,ld,l1,l2,l4)

          la=int((aa(3)-z0)/hdz)
          lb=int((bb(3)-z0)/hdz)
          ld=int((dd(3)-z0)/hdz)
          l1=int((f1(3)-z0)/hdz)
          l2=int((f2(3)-z0)/hdz)
          l4=int((f4(3)-z0)/hdz)
          ku=max0(la,lb,ld,l1,l2,l4)
          kd=min0(la,lb,ld,l1,l2,l4)

          xn1=-(f2(2)-f1(2))*(f4(3)-f1(3))+(f2(3)-f1(3))*(f4(2)-f1(2))
          yn1=-(f2(3)-f1(3))*(f4(1)-f1(1))+(f2(1)-f1(1))*(f4(3)-f1(3))
          zn1=-(f2(1)-f1(1))*(f4(2)-f1(2))+(f2(2)-f1(2))*(f4(1)-f1(1))

          xn2=(bb(2)-aa(2))*(dd(3)-aa(3))-(bb(3)-aa(3))*(dd(2)-aa(2))
          yn2=(bb(3)-aa(3))*(dd(1)-aa(1))-(bb(1)-aa(1))*(dd(3)-aa(3))
          zn2=(bb(1)-aa(1))*(dd(2)-aa(2))-(bb(2)-aa(2))*(dd(1)-aa(1))

          xn3=-(f4(2)-f1(2))*(aa(3)-f1(3))+(f4(3)-f1(3))*(aa(2)-f1(2))
          yn3=-(f4(3)-f1(3))*(aa(1)-f1(1))+(f4(1)-f1(1))*(aa(3)-f1(3))
          zn3=-(f4(1)-f1(1))*(aa(2)-f1(2))+(f4(2)-f1(2))*(aa(1)-f1(1))

          xn4=-(aa(2)-f1(2))*(f2(3)-f1(3))+(aa(3)-f1(3))*(f2(2)-f1(2))
          yn4=-(aa(3)-f1(3))*(f2(1)-f1(1))+(aa(1)-f1(1))*(f2(3)-f1(3))
          zn4=-(aa(1)-f1(1))*(f2(2)-f1(2))+(aa(2)-f1(2))*(f2(1)-f1(1))

          xn5=-(bb(2)-f2(2))*(f4(3)-f2(3))+(bb(3)-f2(3))*(f4(2)-f2(2))
          yn5=-(bb(3)-f2(3))*(f4(1)-f2(1))+(bb(1)-f2(1))*(f4(3)-f2(3))
          zn5=-(bb(1)-f2(1))*(f4(2)-f2(2))+(bb(2)-f2(2))*(f4(1)-f2(1))

          do k=kd,ku
            do j=jd,ju
              do i=id,iu                
                side1=xn1*(x(i)-f1(1))+yn1*(y(j)-f1(2))+zn1*(z(k)-f1(3))
                side2=xn2*(x(i)-aa(1))+yn2*(y(j)-aa(2))+zn2*(z(k)-aa(3))
                side3=xn3*(x(i)-f1(1))+yn3*(y(j)-f1(2))+zn3*(z(k)-f1(3))
                side4=xn4*(x(i)-f1(1))+yn4*(y(j)-f1(2))+zn4*(z(k)-f1(3))
                side5=xn5*(x(i)-f2(1))+yn5*(y(j)-f2(2))+zn5*(z(k)-f2(3))
                if(side1.ge.0.0d0.and.side2.ge.0.0d0.and.side3.ge.0.0d0
     .          .and.side4.ge.0.0d0.and.side5.ge.0.0d0) then
                  if(mod(i,2).eq.0) then
                    ic=i/2+1
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).ne.0) then
                        ke=(k+1)/2
                        w(ic,jc,ke)=zsct(m)+omegax(m)*(y(j)-ysc(m))-
     .                                      omegay(m)*(x(i)-xsc(m))
                      endif
                    else
                      je=(j+1)/2
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        v(ic,je,kc)=ysct(m)+omegaz(m)*(x(i)-xsc(m))-
     .                                      omegax(m)*(z(k)-zsc(m))
                      endif
                    endif
                  else
                    ie=(i+1)/2
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        u(ie,jc,kc)=xsct(m)+omegay(m)*(z(k)-zsc(m))-
     .                                      omegaz(m)*(y(j)-ysc(m))
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
          enddo

c  C-3, B-2, D-4:

          lc=int((cc(1)-x0)/hdx)
          lb=int((bb(1)-x0)/hdx)
          ld=int((dd(1)-x0)/hdx)
          l3=int((f3(1)-x0)/hdx)
          l2=int((f2(1)-x0)/hdx)
          l4=int((f4(1)-x0)/hdx)
          iu=max0(lc,lb,ld,l3,l2,l4)
          id=min0(lc,lb,ld,l3,l2,l4)

          lc=int((cc(2)-y0)/hdy)
          lb=int((bb(2)-y0)/hdy)
          ld=int((dd(2)-y0)/hdy)
          l3=int((f3(2)-y0)/hdy)
          l2=int((f2(2)-y0)/hdy)
          l4=int((f4(2)-y0)/hdy)
          ju=max0(lc,lb,ld,l3,l2,l4)
          jd=min0(lc,lb,ld,l3,l2,l4)

          lc=int((cc(3)-z0)/hdz)
          lb=int((bb(3)-z0)/hdz)
          ld=int((dd(3)-z0)/hdz)
          l3=int((f3(3)-z0)/hdz)
          l2=int((f2(3)-z0)/hdz)
          l4=int((f4(3)-z0)/hdz)
          ku=max0(lc,lb,ld,l3,l2,l4)
          kd=min0(lc,lb,ld,l3,l2,l4)

          xn1=-(f3(2)-f4(2))*(f3(3)-f2(3))+(f3(3)-f4(3))*(f3(2)-f2(2))
          yn1=-(f3(3)-f4(3))*(f3(1)-f2(1))+(f3(1)-f4(1))*(f3(3)-f2(3))
          zn1=-(f3(1)-f4(1))*(f3(2)-f2(2))+(f3(2)-f4(2))*(f3(1)-f2(1))

          xn2=(cc(2)-dd(2))*(cc(3)-bb(3))-(cc(3)-dd(3))*(cc(2)-bb(2))
          yn2=(cc(3)-dd(3))*(cc(1)-bb(1))-(cc(1)-dd(1))*(cc(3)-bb(3))
          zn2=(cc(1)-dd(1))*(cc(2)-bb(2))-(cc(2)-dd(2))*(cc(1)-bb(1))

          xn3=-(f3(2)-cc(2))*(f3(3)-f4(3))+(f3(3)-cc(3))*(f3(2)-f4(2))
          yn3=-(f3(3)-cc(3))*(f3(1)-f4(1))+(f3(1)-cc(1))*(f3(3)-f4(3))
          zn3=-(f3(1)-cc(1))*(f3(2)-f4(2))+(f3(2)-cc(2))*(f3(1)-f4(1))

          xn4=-(f3(2)-f2(2))*(f3(3)-cc(3))+(f3(3)-f2(3))*(f3(2)-cc(2))
          yn4=-(f3(3)-f2(3))*(f3(1)-cc(1))+(f3(1)-f2(1))*(f3(3)-cc(3))
          zn4=-(f3(1)-f2(1))*(f3(2)-cc(2))+(f3(2)-f2(2))*(f3(1)-cc(1))

          xn5=-(f4(2)-f2(2))*(bb(3)-f2(3))+(f4(3)-f2(3))*(bb(2)-f2(2))
          yn5=-(f4(3)-f2(3))*(bb(1)-f2(1))+(f4(1)-f2(1))*(bb(3)-f2(3))
          zn5=-(f4(1)-f2(1))*(bb(2)-f2(2))+(f4(2)-f2(2))*(bb(1)-f2(1))

          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                side1=xn1*(x(i)-f3(1))+yn1*(y(j)-f3(2))+zn1*(z(k)-f3(3))
                side2=xn2*(x(i)-cc(1))+yn2*(y(j)-cc(2))+zn2*(z(k)-cc(3))
                side3=xn3*(x(i)-f3(1))+yn3*(y(j)-f3(2))+zn3*(z(k)-f3(3))
                side4=xn4*(x(i)-f3(1))+yn4*(y(j)-f3(2))+zn4*(z(k)-f3(3))
                side5=xn5*(x(i)-f2(1))+yn5*(y(j)-f2(2))+zn5*(z(k)-f2(3))
                if(side1.ge.0.0d0.and.side2.ge.0.0d0.and.side3.ge.0.0d0
     .          .and.side4.ge.0.0d0.and.side5.ge.0.0d0) then
                  if(mod(i,2).eq.0) then
                    ic=i/2+1
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).ne.0) then
                        ke=(k+1)/2
                        w(ic,jc,ke)=zsct(m)+omegax(m)*(y(j)-ysc(m))-
     .                                      omegay(m)*(x(i)-xsc(m))
                      endif
                    else
                      je=(j+1)/2
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        v(ic,je,kc)=ysct(m)+omegaz(m)*(x(i)-xsc(m))-
     .                                      omegax(m)*(z(k)-zsc(m))
                      endif
                    endif
                  else
                    ie=(i+1)/2
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        u(ie,jc,kc)=xsct(m)+omegay(m)*(z(k)-zsc(m))-
     .                                      omegaz(m)*(y(j)-ysc(m))
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
          enddo

        enddo
      enddo

      IF(i01(m).eq.1) THEN

        do m2=1,ns2-2

c  North pole

          f1(1)=xs(0,0,m)
          f2(1)=xs(0,m2,m)
          f4(1)=xs(0,m2+1,m)

          f1(2)=ys(0,0,m)
          f2(2)=ys(0,m2,m)
          f4(2)=ys(0,m2+1,m)

          f1(3)=zs(0,0,m)
          f2(3)=zs(0,m2,m)
          f4(3)=zs(0,m2+1,m)

          gacobi=fst(0,0,10,m)
          xn=-fst(0,0,7,m)/gacobi
          yn=-fst(0,0,8,m)/gacobi
          zn=-fst(0,0,9,m)/gacobi
          aa(1)=f1(1)+sn*xn
          aa(2)=f1(2)+sn*yn
          aa(3)=f1(3)+sn*zn

          gacobi=fst(0,m2,10,m)
          xn=-fst(0,m2,7,m)/gacobi
          yn=-fst(0,m2,8,m)/gacobi
          zn=-fst(0,m2,9,m)/gacobi
          bb(1)=f2(1)+sn*xn
          bb(2)=f2(2)+sn*yn
          bb(3)=f2(3)+sn*zn

          gacobi=fst(0,m2+1,10,m)
          xn=-fst(0,m2+1,7,m)/gacobi
          yn=-fst(0,m2+1,8,m)/gacobi
          zn=-fst(0,m2+1,9,m)/gacobi
          dd(1)=f4(1)+sn*xn
          dd(2)=f4(2)+sn*yn
          dd(3)=f4(3)+sn*zn

          la=int((aa(1)-x0)/hdx)
          lb=int((bb(1)-x0)/hdx)
          ld=int((dd(1)-x0)/hdx)
          l1=int((f1(1)-x0)/hdx)
          l2=int((f2(1)-x0)/hdx)
          l4=int((f4(1)-x0)/hdx)
          iu=max0(la,lb,ld,l1,l2,l4)
          id=min0(la,lb,ld,l1,l2,l4)

          la=int((aa(2)-y0)/hdy)
          lb=int((bb(2)-y0)/hdy)
          ld=int((dd(2)-y0)/hdy)
          l1=int((f1(2)-y0)/hdy)
          l2=int((f2(2)-y0)/hdy)
          l4=int((f4(2)-y0)/hdy)
          ju=max0(la,lb,ld,l1,l2,l4)
          jd=min0(la,lb,ld,l1,l2,l4)

          la=int((aa(3)-z0)/hdz)
          lb=int((bb(3)-z0)/hdz)
          ld=int((dd(3)-z0)/hdz)
          l1=int((f1(3)-z0)/hdz)
          l2=int((f2(3)-z0)/hdz)
          l4=int((f4(3)-z0)/hdz)
          ku=max0(la,lb,ld,l1,l2,l4)
          kd=min0(la,lb,ld,l1,l2,l4)

          xn1=-(f2(2)-f1(2))*(f4(3)-f1(3))+(f2(3)-f1(3))*(f4(2)-f1(2))
          yn1=-(f2(3)-f1(3))*(f4(1)-f1(1))+(f2(1)-f1(1))*(f4(3)-f1(3))
          zn1=-(f2(1)-f1(1))*(f4(2)-f1(2))+(f2(2)-f1(2))*(f4(1)-f1(1))

          xn2=(bb(2)-aa(2))*(dd(3)-aa(3))-(bb(3)-aa(3))*(dd(2)-aa(2))
          yn2=(bb(3)-aa(3))*(dd(1)-aa(1))-(bb(1)-aa(1))*(dd(3)-aa(3))
          zn2=(bb(1)-aa(1))*(dd(2)-aa(2))-(bb(2)-aa(2))*(dd(1)-aa(1))

          xn3=-(f4(2)-f1(2))*(aa(3)-f1(3))+(f4(3)-f1(3))*(aa(2)-f1(2))
          yn3=-(f4(3)-f1(3))*(aa(1)-f1(1))+(f4(1)-f1(1))*(aa(3)-f1(3))
          zn3=-(f4(1)-f1(1))*(aa(2)-f1(2))+(f4(2)-f1(2))*(aa(1)-f1(1))

          xn4=-(aa(2)-f1(2))*(f2(3)-f1(3))+(aa(3)-f1(3))*(f2(2)-f1(2))
          yn4=-(aa(3)-f1(3))*(f2(1)-f1(1))+(aa(1)-f1(1))*(f2(3)-f1(3))
          zn4=-(aa(1)-f1(1))*(f2(2)-f1(2))+(aa(2)-f1(2))*(f2(1)-f1(1))

          xn5=-(bb(2)-f2(2))*(f4(3)-f2(3))+(bb(3)-f2(3))*(f4(2)-f2(2))
          yn5=-(bb(3)-f2(3))*(f4(1)-f2(1))+(bb(1)-f2(1))*(f4(3)-f2(3))
          zn5=-(bb(1)-f2(1))*(f4(2)-f2(2))+(bb(2)-f2(2))*(f4(1)-f2(1))

          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                side1=xn1*(x(i)-f1(1))+yn1*(y(j)-f1(2))+zn1*(z(k)-f1(3))
                side2=xn2*(x(i)-aa(1))+yn2*(y(j)-aa(2))+zn2*(z(k)-aa(3))
                side3=xn3*(x(i)-f1(1))+yn3*(y(j)-f1(2))+zn3*(z(k)-f1(3))
                side4=xn4*(x(i)-f1(1))+yn4*(y(j)-f1(2))+zn4*(z(k)-f1(3))
                side5=xn5*(x(i)-f2(1))+yn5*(y(j)-f2(2))+zn5*(z(k)-f2(3))
                if(side1.ge.0.0d0.and.side2.ge.0.0d0.and.side3.ge.0.0d0
     .          .and.side4.ge.0.0d0.and.side5.ge.0.0d0) then
                  if(mod(i,2).eq.0) then
                    ic=i/2+1
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).ne.0) then
                        ke=(k+1)/2
                        w(ic,jc,ke)=zsct(m)+omegax(m)*(y(j)-ysc(m))-
     .                                      omegay(m)*(x(i)-xsc(m))
                      endif
                    else
                      je=(j+1)/2
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        v(ic,je,kc)=ysct(m)+omegaz(m)*(x(i)-xsc(m))-
     .                                      omegax(m)*(z(k)-zsc(m))
                      endif
                    endif
                  else
                    ie=(i+1)/2
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        u(ie,jc,kc)=xsct(m)+omegay(m)*(z(k)-zsc(m))-
     .                                      omegaz(m)*(y(j)-ysc(m))
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
          enddo

c  South pole

          f3(1)=xs(ns1,0,m)
          f4(1)=xs(ns1,m2+1,m)
          f2(1)=xs(ns1,m2,m)

          f3(2)=ys(ns1,0,m)
          f4(2)=ys(ns1,m2+1,m)
          f2(2)=ys(ns1,m2,m)

          f3(3)=zs(ns1,0,m)
          f4(3)=zs(ns1,m2+1,m)
          f2(3)=zs(ns1,m2,m)

          gacobi=fst(ns1,0,10,m)
          xn=-fst(ns1,0,7,m)/gacobi
          yn=-fst(ns1,0,8,m)/gacobi
          zn=-fst(ns1,0,9,m)/gacobi
          cc(1)=f3(1)+sn*xn
          cc(2)=f3(2)+sn*yn
          cc(3)=f3(3)+sn*zn

          gacobi=fst(ns1,m2+1,10,m)
          xn=-fst(ns1,m2+1,7,m)/gacobi
          yn=-fst(ns1,m2+1,8,m)/gacobi
          zn=-fst(ns1,m2+1,9,m)/gacobi
          dd(1)=f4(1)+sn*xn
          dd(2)=f4(2)+sn*yn
          dd(3)=f4(3)+sn*zn

          gacobi=fst(ns1,m2,10,m)
          xn=-fst(ns1,m2,7,m)/gacobi
          yn=-fst(ns1,m2,8,m)/gacobi
          zn=-fst(ns1,m2,9,m)/gacobi
          bb(1)=f2(1)+sn*xn
          bb(2)=f2(2)+sn*yn
          bb(3)=f2(3)+sn*zn

          lc=int((cc(1)-x0)/hdx)
          ld=int((dd(1)-x0)/hdx)
          lb=int((bb(1)-x0)/hdx)
          l3=int((f3(1)-x0)/hdx)
          l4=int((f4(1)-x0)/hdx)
          l2=int((f2(1)-x0)/hdx)
          iu=max0(lc,ld,lb,l3,l4,l2)
          id=min0(lc,ld,lb,l3,l4,l2)

          lc=int((cc(2)-y0)/hdy)
          ld=int((dd(2)-y0)/hdy)
          lb=int((bb(2)-y0)/hdy)
          l3=int((f3(2)-y0)/hdy)
          l4=int((f4(2)-y0)/hdy)
          l2=int((f2(2)-y0)/hdy)
          ju=max0(lc,ld,lb,l3,l4,l2)
          jd=min0(lc,ld,lb,l3,l4,l2)

          lc=int((cc(3)-z0)/hdz)
          ld=int((dd(3)-z0)/hdz)
          lb=int((bb(3)-z0)/hdz)
          l3=int((f3(3)-z0)/hdz)
          l4=int((f4(3)-z0)/hdz)
          l2=int((f2(3)-z0)/hdz)
          ku=max0(lc,ld,lb,l3,l4,l2)
          kd=min0(lc,ld,lb,l3,l4,l2)

          xn1=-(f3(2)-f4(2))*(f3(3)-f2(3))+(f3(3)-f4(3))*(f3(2)-f2(2))
          yn1=-(f3(3)-f4(3))*(f3(1)-f2(1))+(f3(1)-f4(1))*(f3(3)-f2(3))
          zn1=-(f3(1)-f4(1))*(f3(2)-f2(2))+(f3(2)-f4(2))*(f3(1)-f2(1))

          xn2=(cc(2)-dd(2))*(cc(3)-bb(3))-(cc(3)-dd(3))*(cc(2)-bb(2))
          yn2=(cc(3)-dd(3))*(cc(1)-bb(1))-(cc(1)-dd(1))*(cc(3)-bb(3))
          zn2=(cc(1)-dd(1))*(cc(2)-bb(2))-(cc(2)-dd(2))*(cc(1)-bb(1))

          xn3=-(f3(2)-cc(2))*(f3(3)-f4(3))+(f3(3)-cc(3))*(f3(2)-f4(2))
          yn3=-(f3(3)-cc(3))*(f3(1)-f4(1))+(f3(1)-cc(1))*(f3(3)-f4(3))
          zn3=-(f3(1)-cc(1))*(f3(2)-f4(2))+(f3(2)-cc(2))*(f3(1)-f4(1))

          xn4=-(f3(2)-f2(2))*(f3(3)-cc(3))+(f3(3)-f2(3))*(f3(2)-cc(2))
          yn4=-(f3(3)-f2(3))*(f3(1)-cc(1))+(f3(1)-f2(1))*(f3(3)-cc(3))
          zn4=-(f3(1)-f2(1))*(f3(2)-cc(2))+(f3(2)-f2(2))*(f3(1)-cc(1))

          xn5=-(f4(2)-f2(2))*(bb(3)-f2(3))+(f4(3)-f2(3))*(bb(2)-f2(2))
          yn5=-(f4(3)-f2(3))*(bb(1)-f2(1))+(f4(1)-f2(1))*(bb(3)-f2(3))
          zn5=-(f4(1)-f2(1))*(bb(2)-f2(2))+(f4(2)-f2(2))*(bb(1)-f2(1))

          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                side1=xn1*(x(i)-f3(1))+yn1*(y(j)-f3(2))+zn1*(z(k)-f3(3))
                side2=xn2*(x(i)-cc(1))+yn2*(y(j)-cc(2))+zn2*(z(k)-cc(3))
                side3=xn3*(x(i)-f3(1))+yn3*(y(j)-f3(2))+zn3*(z(k)-f3(3))
                side4=xn4*(x(i)-f3(1))+yn4*(y(j)-f3(2))+zn4*(z(k)-f3(3))
                side5=xn5*(x(i)-f2(1))+yn5*(y(j)-f2(2))+zn5*(z(k)-f2(3))
                if(side1.ge.0.0d0.and.side2.ge.0.0d0.and.side3.ge.0.0d0
     .          .and.side4.ge.0.0d0.and.side5.ge.0.0d0) then
                  if(mod(i,2).eq.0) then
                    ic=i/2+1
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).ne.0) then
                        ke=(k+1)/2
                        w(ic,jc,ke)=zsct(m)+omegax(m)*(y(j)-ysc(m))-
     .                                      omegay(m)*(x(i)-xsc(m))
                      endif
                    else
                      je=(j+1)/2
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        v(ic,je,kc)=ysct(m)+omegaz(m)*(x(i)-xsc(m))-
     .                                      omegax(m)*(z(k)-zsc(m))
                      endif
                    endif
                  else
                    ie=(i+1)/2
                    if(mod(j,2).eq.0) then
                      jc=j/2+1
                      if(mod(k,2).eq.0) then
                        kc=k/2+1
                        u(ie,jc,kc)=xsct(m)+omegay(m)*(z(k)-zsc(m))-
     .                                      omegaz(m)*(y(j)-ysc(m))
                      endif
                    endif
                  endif
                endif
              enddo
            enddo
          enddo

        enddo

      ENDIF

      ENDDO

      return
      end


c-----------------------------------------------------------------------
