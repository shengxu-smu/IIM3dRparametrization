c-----------------------------------------------------------------------
c
      subroutine singular_force
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer isphr,nijk,nn,id,jd,kd,ic,ie,jc,je,kc,ke,ifilter
      integer iu,ju,ku,iv,jv,kv,iw,jw,kw
      parameter(nijk=3,nn=3)
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,xe1,ye1,ze1,xe2,ye2,ze2
      real*8 sn,xx,yy,zz,sx,sy,sz,su,sv,sw,dudn,dvdn,dwdn
      real*8 foo,os1,os2,os3,ost1,ost2,osnn
      real*8 gacobi,g1,g2,g12,xn1,yn1,zn1,xn2,yn2,zn2,dot1dn,dot2dn
      real*8 hi(nijk),fi(nijk),hj(nijk),fj(nijk),hk(nijk),fk(nijk)
      real*8 h1(nijk),f1(nijk),h2(nijk),f2(nijk),h3(nijk),f3(nijk)
      real*8 uu(nn),vv(nn),ww(nn),o1(nn),o2(nn),o3(nn),ot1(nn),ot2(nn)
      real*8 rhs1(0:ns1,0:ns2),rhs2(0:ns1,0:ns2),rhs0(0:ns1,0:ns2)
      real*8 fs(2*ns1+2,ns2),ft(ns1,ns2)

      ifilter=0

      do k=0,nz
        do j=0,ny
          do i=0,nx
            fcee(i,j,k)=(w(i,j+1,k)-w(i,j,k))*dy1-
     .                  (v(i,j,k+1)-v(i,j,k))*dz1
            fece(i,j,k)=(u(i,j,k+1)-u(i,j,k))*dz1-
     .                  (w(i+1,j,k)-w(i,j,k))*dx1
            feec(i,j,k)=(v(i+1,j,k)-v(i,j,k))*dx1-
     .                  (u(i,j+1,k)-u(i,j,k))*dy1
          enddo
        enddo
      enddo

      sn=1.01d0*sqrt(dx*dx+dy*dy+dz*dz)

      DO m=1,ms

      isphr=i01(m)

      do m2=0,ns2-1
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
          g1=xt*xt+yt*yt+zt*zt
          g2=xb*xb+yb*yb+zb*zb
          g12=xt*xb+yt*yb+zt*zb
          gacobi=fst(m1,m2,10,m)
          xn=xn/gacobi
          yn=yn/gacobi
          zn=zn/gacobi
          xe1=yb*zn-zb*yn
          ye1=zb*xn-xb*zn
          ze1=xb*yn-yb*xn
          xe2=yn*zt-zn*yt
          ye2=zn*xt-xn*zt
          ze2=xn*yt-yn*xt

          id=int(sign(1.0,xn))
          jd=int(sign(1.0,yn))
          kd=int(sign(1.0,zn))

          sx=xs(m1,m2,m)
          sy=ys(m1,m2,m)
          sz=zs(m1,m2,m)
          su=us(m1,m2,m)
          sv=vs(m1,m2,m)
          sw=ws(m1,m2,m)

          xn1=snd(m1,m2,10,m)
          yn1=snd(m1,m2,11,m)
          zn1=snd(m1,m2,12,m)
          xn2=snd(m1,m2,13,m)
          yn2=snd(m1,m2,14,m)
          zn2=snd(m1,m2,15,m)

          do n=1,nn

            xx=sx+sn*xn*dble(n)
            yy=sy+sn*yn*dble(n)
            zz=sz+sn*zn*dble(n)
            i=int((xx-x0)/hdx)
            j=int((yy-y0)/hdy)
            k=int((zz-z0)/hdz)
            if(mod(i,2).eq.0) then
              ic=i/2+1
              ie=ic-1
            else
              ie=(i+1)/2
              ic=ie
            endif
            if(mod(j,2).eq.0) then
              jc=j/2+1
              je=jc-1
            else
              je=(j+1)/2
              jc=je
            endif
            if(mod(k,2).eq.0) then
              kc=k/2+1
              ke=kc-1
            else
              ke=(k+1)/2
              kc=ke
            endif

            iu=ie
            ju=jc
            ku=kc
            iv=ic
            jv=je
            kv=kc
            iw=ic
            jw=jc
            kw=ke
            i1=ic
            j1=je
            k1=ke
            i2=ie
            j2=jc
            k2=ke
            i3=ie
            j3=je
            k3=kc
            if(id.lt.0.0d0) then
              iu=iu+1
              iv=iv+1
              iw=iw+1
              i1=i1+1
              i2=i2+1
              i3=i3+1
            endif
            if(jd.lt.0.0d0) then
              ju=ju+1
              jv=jv+1
              jw=jw+1
              j1=j1+1
              j2=j2+1
              j3=j3+1
            endif
            if(kd.lt.0.0d0) then
              ku=ku+1
              kv=kv+1
              kw=kw+1
              k1=k1+1
              k2=k2+1
              k3=k3+1
            endif

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xe(iu+id*(i-1))
                  fi(i)=u(iu+id*(i-1),ju+jd*(j-1),ku+kd*(k-1))
                  h1(i)=xc(i1+id*(i-1))
                  f1(i)=fcee(i1+id*(i-1),j1+jd*(j-1),k1+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                call interpolate(h1,f1,nijk,xx,f2(j),foo)
                hj(j)=yc(ju+jd*(j-1))
                h2(j)=ye(j1+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              call interpolate(h2,f2,nijk,yy,f3(k),foo)
              hk(k)=zc(ku+kd*(k-1))
              h3(k)=ze(k1+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,uu(n),foo)
            call interpolate(h3,f3,nijk,zz,o1(n),foo)

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xc(iv+id*(i-1))
                  fi(i)=v(iv+id*(i-1),jv+jd*(j-1),kv+kd*(k-1))
                  h1(i)=xe(i2+id*(i-1))
                  f1(i)=fece(i2+id*(i-1),j2+jd*(j-1),k2+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                call interpolate(h1,f1,nijk,xx,f2(j),foo)
                hj(j)=ye(jv+jd*(j-1))
                h2(j)=yc(j2+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              call interpolate(h2,f2,nijk,yy,f3(k),foo)
              hk(k)=zc(kv+kd*(k-1))
              h3(k)=ze(k2+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,vv(n),foo)
            call interpolate(h3,f3,nijk,zz,o2(n),foo)

            do k=1,nijk
              do j=1,nijk
                do i=1,nijk
                  hi(i)=xc(iw+id*(i-1))
                  fi(i)=w(iw+id*(i-1),jw+jd*(j-1),kw+kd*(k-1))
                  h1(i)=xe(i3+id*(i-1))
                  f1(i)=feec(i3+id*(i-1),j3+jd*(j-1),k3+kd*(k-1))
                enddo
                call interpolate(hi,fi,nijk,xx,fj(j),foo)
                call interpolate(h1,f1,nijk,xx,f2(j),foo)
                hj(j)=yc(jw+jd*(j-1))
                h2(j)=ye(j3+jd*(j-1))
              enddo
              call interpolate(hj,fj,nijk,yy,fk(k),foo)
              call interpolate(h2,f2,nijk,yy,f3(k),foo)
              hk(k)=ze(kw+kd*(k-1))
              h3(k)=zc(k3+kd*(k-1))
            enddo
            call interpolate(hk,fk,nijk,zz,ww(n),foo)
            call interpolate(h3,f3,nijk,zz,o3(n),foo)

            ot1(n)=o1(n)*xt+o2(n)*yt+o3(n)*zt
            ot2(n)=o1(n)*xb+o2(n)*yb+o3(n)*zb

          enddo

          if(nn.eq.1) then
            dudn=(uu(1)-su)/sn
            dvdn=(vv(1)-sv)/sn
            dwdn=(ww(1)-sw)/sn
          endif
          if(nn.eq.2) then
            dudn=(4.0d0*uu(1)-uu(2)-3.0d0*su)/(2.0d0*sn)
            dvdn=(4.0d0*vv(1)-vv(2)-3.0d0*sv)/(2.0d0*sn)
            dwdn=(4.0d0*ww(1)-ww(2)-3.0d0*sw)/(2.0d0*sn)
          endif
          if(nn.eq.3) then
            dudn=(18.0d0*uu(1)-9.0d0*uu(2)+2.0d0*uu(3)-11.0d0*su)
     .           /(6.0d0*sn)
            dvdn=(18.0d0*vv(1)-9.0d0*vv(2)+2.0d0*vv(3)-11.0d0*sv)
     .           /(6.0d0*sn)
            dwdn=(18.0d0*ww(1)-9.0d0*ww(2)+2.0d0*ww(3)-11.0d0*sw)
     .           /(6.0d0*sn)
          endif

          fr(m1,m2,1,m)=gacobi*(dudn-(omegay(m)*zn-omegaz(m)*yn))
          fr(m1,m2,2,m)=gacobi*(dvdn-(omegaz(m)*xn-omegax(m)*zn))
          fr(m1,m2,3,m)=gacobi*(dwdn-(omegax(m)*yn-omegay(m)*xn))

          ost1=-(dudn*xe2+dvdn*ye2+dwdn*ze2)+
     .         (omegax(m)*xt+omegay(m)*yt+omegaz(m)*zt)
          ost2=(dudn*xe1+dvdn*ye1+dwdn*ze1)+
     .         (omegax(m)*xb+omegay(m)*yb+omegaz(m)*zb)
          osnn=2.0d0*gacobi*(omegax(m)*xn+omegay(m)*yn+omegaz(m)*zn)

          os1=(xe1*ost1+xe2*ost2+xn*osnn)/gacobi
          os2=(ye1*ost1+ye2*ost2+yn*osnn)/gacobi
          os3=(ze1*ost1+ze2*ost2+zn*osnn)/gacobi

          if(nn.eq.1) then
            dot1dn=(ot1(1)-ost1)/sn
            dot2dn=(ot2(1)-ost2)/sn
          endif
          if(nn.eq.2) then
            dot1dn=(4.0d0*ot1(1)-ot1(2)-3.0d0*ost1)/(2.0d0*sn)
            dot2dn=(4.0d0*ot2(1)-ot2(2)-3.0d0*ost2)/(2.0d0*sn)
          endif
          if(nn.eq.3) then
            dot1dn=(18.0d0*ot1(1)-9.0d0*ot1(2)+2.0d0*ot1(3)-11.0d0*ost1)
     .           /(6.0d0*sn)
            dot2dn=(18.0d0*ot2(1)-9.0d0*ot2(2)+2.0d0*ot2(3)-11.0d0*ost2)
     .           /(6.0d0*sn)
          endif

          rhs1(m1,m2)=-(((2.0d0*omegax(m)-os1)*(g1*xn2-g12*xn1)+
     .                   (2.0d0*omegay(m)-os2)*(g1*yn2-g12*yn1)+
     .                   (2.0d0*omegaz(m)-os3)*(g1*zn2-g12*zn1))/gacobi+
     .                  g12*dot1dn-g1*dot2dn)/(Re*gacobi)
     .               -xt*(omegayt(m)*(sz-zsc(m))-omegazt(m)*(sy-ysc(m)))
     .               -yt*(omegazt(m)*(sx-xsc(m))-omegaxt(m)*(sz-zsc(m)))
     .               -zt*(omegaxt(m)*(sy-ysc(m))-omegayt(m)*(sx-xsc(m)))
          rhs2(m1,m2)=-(((2.0d0*omegax(m)-os1)*(g12*xn2-g2*xn1)+
     .                   (2.0d0*omegay(m)-os2)*(g12*yn2-g2*yn1)+
     .                   (2.0d0*omegaz(m)-os3)*(g12*zn2-g2*zn1))/gacobi+
     .                  g2*dot1dn-g12*dot2dn)/(Re*gacobi)
     .               -xb*(omegayt(m)*(sz-zsc(m))-omegazt(m)*(sy-ysc(m)))
     .               -yb*(omegazt(m)*(sx-xsc(m))-omegaxt(m)*(sz-zsc(m)))
     .               -zb*(omegaxt(m)*(sy-ysc(m))-omegayt(m)*(sx-xsc(m)))
        enddo
      enddo

      do m1=0,ns1
        fr(m1,ns2,1,m)=fr(m1,0,1,m)
        fr(m1,ns2,2,m)=fr(m1,0,2,m)
        fr(m1,ns2,3,m)=fr(m1,0,3,m)
        rhs1(m1,ns2)=rhs1(m1,0)
        rhs2(m1,ns2)=rhs2(m1,0)
      enddo

      if(isphr.eq.0) then
        do m2=0,ns2
          fr(ns1,m2,1,m)=fr(0,m2,1,m)
          fr(ns1,m2,2,m)=fr(0,m2,2,m)
          fr(ns1,m2,3,m)=fr(0,m2,3,m)
          rhs1(ns1,m2)=rhs1(0,m2)
          rhs2(ns1,m2)=rhs2(0,m2)
        enddo
      endif

      IF(ifilter.eq.1) THEN
      do m2=0,ns2
        do m1=0,ns1
          rhs0(m1,m2)=fr(m1,m2,1,m)
        enddo
      enddo
      call filter(rhs0,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,1,m)=rhs0(m1,m2)
        enddo
      enddo
      do m2=0,ns2
        do m1=0,ns1
          rhs0(m1,m2)=fr(m1,m2,2,m)
        enddo
      enddo
      call filter(rhs0,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,2,m)=rhs0(m1,m2)
        enddo
      enddo
      do m2=0,ns2
        do m1=0,ns1
          rhs0(m1,m2)=fr(m1,m2,3,m)
        enddo
      enddo
      call filter(rhs0,isphr,32,32,1.0d0)
      do m2=0,ns2
        do m1=0,ns1
          fr(m1,m2,3,m)=rhs0(m1,m2)
        enddo
      enddo
      call filter(rhs1,isphr,32,32,-1.0d0)
      call filter(rhs2,isphr,32,32,1.0d0)
      ENDIF

      call spline_derivative(rhs1,rhs0,isphr,1,2)
      call spline_derivative(rhs2,rhs1,isphr,2,2)

      do m2=0,ns2
        do m1=0,ns1
          rhs0(m1,m2)=rhs0(m1,m2)+rhs1(m1,m2)
        enddo
      enddo

      if(isphr.eq.1) then
        do m2=0,ns2-1
          do m1=0,ns1
            fs(m1+1,m2+1)=rhs0(m1,m2)
          enddo
        enddo
        do m2=0,ns2/2
          do m1=ns1+1,2*ns1+1
            fs(m1+1,m2+1)=rhs0(2*ns1+1-m1,ns2/2+m2)
          enddo
        enddo
        do m2=ns2/2+1,ns2-1
          do m1=ns1+1,2*ns1+1
            fs(m1+1,m2+1)=rhs0(2*ns1+1-m1,m2-ns2/2)
          enddo
        enddo

        call poisson2d_sphere(2*ns1+2,ns2,dalfa11,dalfa20,fs)

        do m2=0,ns2-1
          do m1=0,ns1
            fr(m1,m2,4,m)=fs(m1+1,m2+1)
          enddo
        enddo
        do m1=0,ns1
          fr(m1,ns2,4,m)=fr(m1,0,4,m)
        enddo
      else
        do m2=0,ns2-1
          do m1=0,ns1-1
            ft(m1+1,m2+1)=rhs0(m1,m2)
          enddo
        enddo

        call poisson2d_torus(ns1,ns2,dalfa10,dalfa20,ft)

        do m2=0,ns2-1
          do m1=0,ns1-1
            fr(m1,m2,4,m)=ft(m1+1,m2+1)
          enddo
        enddo
        do m1=0,ns1
          fr(m1,ns2,4,m)=fr(m1,0,4,m)
        enddo
        do m2=0,ns2
          fr(ns1,m2,4,m)=fr(0,m2,4,m)
        enddo
      endif

      ENDDO
       
      return
      end


c-----------------------------------------------------------------------
