c-----------------------------------------------------------------------
c
      subroutine fn_correct
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer interp,isphr,nijk,nn,id,jd,kd,ic,jc,kc
      parameter(nijk=2,nn=3)
      real*8 xn,yn,zn,sn,xx,yy,zz,sx,sy,sz,foo,gacobi,dpdn,omega2
      real*8 pplus,pminus
      real*8 h1(nijk),f1(nijk),h2(nijk),f2(nijk),h3(nijk),f3(nijk)
      real*8 pp(nn),qn(0:ns1,0:ns2)

      interp=1
      sn=1.01d0*sqrt(dx*dx+dy*dy+dz*dz)

      IF(interp.eq.1) THEN
      DO m=1,ms

      isphr=i01(m)

      do m2=0,ns2-1
        do m1=0,ns1

          xn=fst(m1,m2,7,m)
          yn=fst(m1,m2,8,m)
          zn=fst(m1,m2,9,m)
          gacobi=fst(m1,m2,10,m)
          xn=xn/gacobi
          yn=yn/gacobi
          zn=zn/gacobi

          sx=xs(m1,m2,m)
          sy=ys(m1,m2,m)
          sz=zs(m1,m2,m)

          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          kd=int(sign(1.0,zn))

          xx=sx+sn*xn
          yy=sy+sn*yn
          zz=sz+sn*zn
          i=int((xx-x0)/hdx)
          j=int((yy-y0)/hdy)
          k=int((zz-z0)/hdz)
          if(mod(i,2).eq.0) then
            ic=i/2+1
          else
            ic=(i+1)/2
          endif
          if(mod(j,2).eq.0) then
            jc=j/2+1
          else
            jc=(j+1)/2
          endif
          if(mod(k,2).eq.0) then
            kc=k/2+1
          else
            kc=(k+1)/2
          endif

          if(id.lt.0.0d0) then
            ic=ic+1
          endif
          if(jd.lt.0.0d0) then
            jc=jc+1
          endif
          if(kd.lt.0.0d0) then
            kc=kc+1
          endif

          do k=1,nijk
            do j=1,nijk
              do i=1,nijk
                h1(i)=xc(ic+id*(i-1))
                f1(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
              enddo
              call interpolate(h1,f1,nijk,xx,f2(j),foo)
              h2(j)=yc(jc+jd*(j-1))
            enddo
            call interpolate(h2,f2,nijk,yy,f3(k),foo)
            h3(k)=zc(kc+kd*(k-1))
          enddo
          call interpolate(h3,f3,nijk,zz,pplus,foo)

          xn=-xn
          yn=-yn
          zn=-zn

          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          kd=int(sign(1.0,zn))

          xx=sx+sn*xn
          yy=sy+sn*yn
          zz=sz+sn*zn
          i=int((xx-x0)/hdx)
          j=int((yy-y0)/hdy)
          k=int((zz-z0)/hdz)

          if(mod(i,2).eq.0) then
            ic=i/2+1
          else
            ic=(i+1)/2
          endif
          if(mod(j,2).eq.0) then
            jc=j/2+1
          else
            jc=(j+1)/2
          endif
          if(mod(k,2).eq.0) then
            kc=k/2+1
          else
            kc=(k+1)/2
          endif

          if(id.lt.0.0d0) then
            ic=ic+1
          endif
          if(jd.lt.0.0d0) then
            jc=jc+1
          endif
          if(kd.lt.0.0d0) then
            kc=kc+1
          endif

          do k=1,nijk
            do j=1,nijk
              do i=1,nijk
                h1(i)=xc(ic+id*(i-1))
                f1(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
              enddo
              call interpolate(h1,f1,nijk,xx,f2(j),foo)
              h2(j)=yc(jc+jd*(j-1))
            enddo
            call interpolate(h2,f2,nijk,yy,f3(k),foo)
            h3(k)=zc(kc+kd*(k-1))
          enddo
          call interpolate(h3,f3,nijk,zz,pminus,foo)

          fr(m1,m2,4,m)=0.5d0*(pplus+pminus+fr(m1,m2,4,m))-
     .                  0.5d0*sn*fr(m1,m2,20,m)/gacobi
        enddo
      enddo

      do m1=0,ns1
        fr(m1,ns2,4,m)=fr(m1,0,4,m)
      enddo
      if(isphr.eq.0) then
        do m2=0,ns2
          fr(ns1,m2,4,m)=fr(0,m2,4,m)
        enddo
      endif

      ENDDO
      ENDIF

cccc

      IF(interp.eq.0) THEN
      DO m=1,ms

      isphr=i01(m)

      do m2=0,ns2-1
        do m1=0,ns1

          xn=fst(m1,m2,7,m)
          yn=fst(m1,m2,8,m)
          zn=fst(m1,m2,9,m)
          gacobi=fst(m1,m2,10,m)
          xn=xn/gacobi
          yn=yn/gacobi
          zn=zn/gacobi

          sx=xs(m1,m2,m)
          sy=ys(m1,m2,m)
          sz=zs(m1,m2,m)

          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          kd=int(sign(1.0,zn))

          do n=1,nn

            xx=sx+sn*xn*dble(n)
            yy=sy+sn*yn*dble(n)
            zz=sz+sn*zn*dble(n)
            i=int((xx-x0)/hdx)
            j=int((yy-y0)/hdy)
            k=int((zz-z0)/hdz)
            if(mod(i,2).eq.0) then
              ic=i/2+1
            else
              ic=(i+1)/2
            endif
            if(mod(j,2).eq.0) then
              jc=j/2+1
            else
              jc=(j+1)/2
            endif
            if(mod(k,2).eq.0) then
              kc=k/2+1
            else
              kc=(k+1)/2
            endif

            if(id.lt.0.0d0) then
              ic=ic+1
            endif
            if(jd.lt.0.0d0) then
              jc=jc+1
            endif
            if(kd.lt.0.0d0) then
              kc=kc+1
            endif

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  h1(i)=xc(ic+id*(i-1))
                  f1(i)=p(ic+id*(i-1),jc+jd*(j-1),kc+kd*(k-1))
                enddo
                call interpolate(h1,f1,nijk,xx,f2(j),foo)
                h2(j)=yc(jc+jd*(j-1))
              enddo
              call interpolate(h2,f2,nijk,yy,f3(k),foo)
              h3(k)=zc(kc+kd*(k-1))
            enddo
            call interpolate(h3,f3,nijk,zz,pp(n),foo)

          enddo

          if(nn.eq.2) then
            fr(m1,m2,4,m)=2.0d0*pp(1)-pp(2)
          endif
          if(nn.eq.3) then
            fr(m1,m2,4,m)=3.0d0*pp(1)-3.0d0*pp(2)+pp(3)
          endif

        enddo
      enddo

      do m1=0,ns1
        fr(m1,ns2,4,m)=fr(m1,0,4,m)
      enddo
      if(isphr.eq.0) then
        do m2=0,ns2
          fr(ns1,m2,4,m)=fr(0,m2,4,m)
        enddo
      endif

      ENDDO
      ENDIF

      return
      end


c----------------------------------------------------------------------
