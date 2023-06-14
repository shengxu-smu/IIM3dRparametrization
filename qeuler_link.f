c-----------------------------------------------------------------------
c
      subroutine qeuler_link
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer l1,l2,l3,iu,id,ju,jd,ku,kd,lc,le
      integer isphr,mm1,mm2,lab,lac,lbc,keep
      integer imod(9),im,iprod
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,xn0,yn0,zn0,sn
      real*8 q1,q2,sx,sy,sz,tx,ty,tz,bx,by,bz,an,bn,cn
      real*8 aa(3),bb(3),cc(3),f1(55),f2(55),f3(55),f4(55),f(55)
      real*8 pnt(2,6)

      DO m=1,ms

      isphr=i01(m)

      niejc(m)=0
      nicje(m)=0
      nicjc(m)=0
      niekc(m)=0
      nicke(m)=0
      nickc(m)=0
      njekc(m)=0
      njcke(m)=0
      njckc(m)=0

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

          do n=1,3
            f1(n+21)=pjc1(m1,m2,n,m)
            f2(n+21)=pjc1(m1+1,m2,n,m)
            f3(n+21)=pjc1(m1+1,m2+1,n,m)
            f4(n+21)=pjc1(m1,m2+1,n,m)
          enddo

          do n=1,6
            f1(n+42)=pjc2(m1,m2,n,m)
            f2(n+42)=pjc2(m1+1,m2,n,m)
            f3(n+42)=pjc2(m1+1,m2+1,n,m)
            f4(n+42)=pjc2(m1,m2+1,n,m)

            f1(n+49)=pjc(m1,m2,n,m)
            f2(n+49)=pjc(m1+1,m2,n,m)
            f3(n+49)=pjc(m1+1,m2+1,n,m)
            f4(n+49)=pjc(m1,m2+1,n,m)
          enddo

          f1(49)=pjc0(m1,m2,m)
          f2(49)=pjc0(m1+1,m2,m)
          f3(49)=pjc0(m1+1,m2+1,m)
          f4(49)=pjc0(m1,m2+1,m)


c  A=1, B=2, C=4:

          aa(1)=f1(1)
          aa(2)=f1(2)
          aa(3)=f1(3)
          bb(1)=f2(1)
          bb(2)=f2(2)
          bb(3)=f2(3)
          cc(1)=f4(1)
          cc(2)=f4(2)
          cc(3)=f4(3)

          l1=int((aa(1)-x0)/hdx)
          l2=int((bb(1)-x0)/hdx)
          l3=int((cc(1)-x0)/hdx)
          call maxmin(l1,l2,l3,iu,id)

          l1=int((aa(2)-y0)/hdy)
          l2=int((bb(2)-y0)/hdy)
          l3=int((cc(2)-y0)/hdy)
          call maxmin(l1,l2,l3,ju,jd)
          
          l1=int((aa(3)-z0)/hdz)
          l2=int((bb(3)-z0)/hdz)
          l3=int((cc(3)-z0)/hdz)
          call maxmin(l1,l2,l3,ku,kd)

          xt=bb(1)-aa(1)
          yt=bb(2)-aa(2)
          zt=bb(3)-aa(3)
          xb=cc(1)-aa(1)
          yb=cc(2)-aa(2)
          zb=cc(3)-aa(3)
          xn=yt*zb-zt*yb
          yn=zt*xb-xt*zb
          zn=xt*yb-yt*xb
          xn0=xn
          yn0=yn
          zn0=zn          
          sn=dsqrt(xn*xn+yn*yn+zn*zn)

c  ij-line

          IF(dabs(zn)/sn.gt.1.0d-12) THEN
          do j=jd,ju
            do i=id,iu
              call line_side(1,2,aa,bb,cc,x(i),y(j),lab)
              call line_side(1,2,aa,cc,bb,x(i),y(j),lac)
              call line_side(1,2,bb,cc,aa,x(i),y(j),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                mm1=m1+1
                mm2=m2-1
                if(m2.eq.0) mm2=ns2-1
                bx=bb(1)-xs(mm1,mm2,m)
                by=bb(2)-ys(mm1,mm2,m)
                cn=xt*by-yt*bx
                if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2+1
                  if(m1.eq.0) mm1=ns1-1
                  tx=cc(1)-xs(mm1,mm2,m)
                  ty=cc(2)-ys(mm1,mm2,m)
                  cn=tx*yb-ty*xb
                  if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2-1
                  if(m1.eq.0) mm1=ns1-1
                  if(m2.eq.0) mm2=ns2-1
                  pnt(1,1)=xt
                  pnt(2,1)=yt
                  pnt(1,2)=xb
                  pnt(2,2)=yb
                  pnt(1,3)=xs(mm1,m2+1,m)-aa(1)
                  pnt(2,3)=ys(mm1,m2+1,m)-aa(2)
                  pnt(1,4)=xs(mm1,m2,m)-aa(1)
                  pnt(2,4)=ys(mm1,m2,m)-aa(2)
                  pnt(1,5)=xs(m1,mm2,m)-aa(1)
                  pnt(2,5)=ys(m1,mm2,m)-aa(2)
                  pnt(1,6)=xs(m1+1,mm2,m)-aa(1)
                  pnt(2,6)=ys(m1+1,mm2,m)-aa(2)
                  call point_polygon(6,pnt,keep)
                endif 
              endif

              if(keep.eq.1) then
                q1=(yb*(x(i)-f1(1))-xb*(y(j)-f1(2)))/zn
                q2=(xt*(y(j)-f1(2))-yt*(x(i)-f1(1)))/zn
                sz=(f2(3)-f1(3))*q1+(f4(3)-f1(3))*q2+f1(3)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif 
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sz
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f2(n)-f1(n))*q1+(f4(n)-f1(n))*q2+f1(n)
                enddo

                k=int((sz-z0)/hdz)
                if(sz.eq.z(k).and.sign(1.0,zn).gt.0.0) k=k-1
                if(mod(k,2).eq.0) then
                  lc=k/2+1
                  le=lc-1
                else
                  le=(k+1)/2
                  lc=le
                endif
                if(mod(i,2).ne.0) then
                  if(mod(j,2).eq.0) then
                    iejc(1,niejc(m),m)=(i+1)/2
                    iejc(2,niejc(m),m)=j/2+1
                    iejc(3,niejc(m),m)=le
                    iejc(4,niejc(m),m)=lc
                    do n=1,55
                      fiejc(n,niejc(m),m)=f(n)
                    enddo
                    niejc(m)=niejc(m)+1
                  endif
                else
                  if(mod(j,2).ne.0) then
                    icje(1,nicje(m),m)=i/2+1
                    icje(2,nicje(m),m)=(j+1)/2
                    icje(3,nicje(m),m)=le
                    icje(4,nicje(m),m)=lc
                    do n=1,55
                      ficje(n,nicje(m),m)=f(n)
                    enddo
                    nicje(m)=nicje(m)+1
                  else
                    icjc(1,nicjc(m),m)=i/2+1
                    icjc(2,nicjc(m),m)=j/2+1
                    icjc(3,nicjc(m),m)=le
                    icjc(4,nicjc(m),m)=lc
                    do n=1,55
                      ficjc(n,nicjc(m),m)=f(n)
                    enddo
                    nicjc(m)=nicjc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

c  ik-line

          IF(dabs(yn)/sn.gt.1.0d-12) THEN
          do k=kd,ku
            do i=id,iu
              call line_side(1,3,aa,bb,cc,x(i),z(k),lab)
              call line_side(1,3,aa,cc,bb,x(i),z(k),lac)
              call line_side(1,3,bb,cc,aa,x(i),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                mm1=m1+1
                mm2=m2-1
                if(m2.eq.0) mm2=ns2-1
                bx=f2(1)-xs(mm1,mm2,m)
                bz=f2(3)-zs(mm1,mm2,m)
                bn=zt*bx-xt*bz
                if(bn*(sign(1.0,yn)+sign(1.0,bn)).ne.0.0) keep=1
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2+1
                  if(m1.eq.0) mm1=ns1-1
                  tx=f4(1)-xs(mm1,mm2,m)
                  tz=f4(3)-zs(mm1,mm2,m)
                  bn=tz*xb-tx*zb
                  if(bn*(sign(1.0,yn)+sign(1.0,bn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2-1
                  if(m1.eq.0) mm1=ns1-1
                  if(m2.eq.0) mm2=ns2-1
                  pnt(1,1)=xt
                  pnt(2,1)=zt
                  pnt(1,2)=xb
                  pnt(2,2)=zb
                  pnt(1,3)=xs(mm1,m2+1,m)-aa(1)
                  pnt(2,3)=zs(mm1,m2+1,m)-aa(3)
                  pnt(1,4)=xs(mm1,m2,m)-aa(1)
                  pnt(2,4)=zs(mm1,m2,m)-aa(3)
                  pnt(1,5)=xs(m1,mm2,m)-aa(1)
                  pnt(2,5)=zs(m1,mm2,m)-aa(3)
                  pnt(1,6)=xs(m1+1,mm2,m)-aa(1)
                  pnt(2,6)=zs(m1+1,mm2,m)-aa(3)
                  call point_polygon(6,pnt,keep)
                endif
              endif

              if(keep.eq.1) then
                q1=-(zb*(x(i)-f1(1))-xb*(z(k)-f1(3)))/yn
                q2=-(xt*(z(k)-f1(3))-zt*(x(i)-f1(1)))/yn
                sy=(f2(2)-f1(2))*q1+(f4(2)-f1(2))*q2+f1(2)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sy
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f2(n)-f1(n))*q1+(f4(n)-f1(n))*q2+f1(n)
                enddo

                j=int((sy-y0)/hdy)
                if(sy.eq.y(j).and.sign(1.0,yn).gt.0.0) j=j-1
                if(mod(j,2).eq.0) then
                  lc=j/2+1
                  le=lc-1
                else
                  le=(j+1)/2
                  lc=le
                endif
                if(mod(i,2).ne.0) then
                  if(mod(k,2).eq.0) then
                    iekc(1,niekc(m),m)=(i+1)/2
                    iekc(2,niekc(m),m)=k/2+1
                    iekc(3,niekc(m),m)=le
                    iekc(4,niekc(m),m)=lc
                    do n=1,55
                      fiekc(n,niekc(m),m)=f(n)
                    enddo
                    niekc(m)=niekc(m)+1
                  endif
                else
                  if(mod(k,2).ne.0) then
                    icke(1,nicke(m),m)=i/2+1
                    icke(2,nicke(m),m)=(k+1)/2
                    icke(3,nicke(m),m)=le
                    icke(4,nicke(m),m)=lc
                    do n=1,55
                      ficke(n,nicke(m),m)=f(n)
                    enddo
                    nicke(m)=nicke(m)+1
                  else
                    ickc(1,nickc(m),m)=i/2+1
                    ickc(2,nickc(m),m)=k/2+1
                    ickc(3,nickc(m),m)=le
                    ickc(4,nickc(m),m)=lc
                    do n=1,55
                      fickc(n,nickc(m),m)=f(n)
                    enddo
                    nickc(m)=nickc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

c  jk-line

          IF(dabs(xn)/sn.gt.1.0d-12) THEN
          do k=kd,ku
            do j=jd,ju
              call line_side(2,3,aa,bb,cc,y(j),z(k),lab)
              call line_side(2,3,aa,cc,bb,y(j),z(k),lac)
              call line_side(2,3,bb,cc,aa,y(j),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                mm1=m1+1
                mm2=m2-1
                if(m2.eq.0) mm2=ns2-1
                by=f2(2)-ys(mm1,mm2,m)
                bz=f2(3)-zs(mm1,mm2,m)
                an=yt*bz-zt*by
                if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2+1
                  if(m1.eq.0) mm1=ns1-1
                  ty=f4(2)-ys(mm1,mm2,m)
                  tz=f4(3)-zs(mm1,mm2,m)
                  an=ty*zb-tz*yb
                  if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(isphr.eq.1.and.m1.eq.0) then
                  keep=0
                else
                  mm1=m1-1
                  mm2=m2-1
                  if(m1.eq.0) mm1=ns1-1
                  if(m2.eq.0) mm2=ns2-1
                  pnt(1,1)=yt
                  pnt(2,1)=zt
                  pnt(1,2)=yb
                  pnt(2,2)=zb
                  pnt(1,3)=ys(mm1,m2+1,m)-aa(2)
                  pnt(2,3)=zs(mm1,m2+1,m)-aa(3)
                  pnt(1,4)=ys(mm1,m2,m)-aa(2)
                  pnt(2,4)=zs(mm1,m2,m)-aa(3)
                  pnt(1,5)=ys(m1,mm2,m)-aa(2)
                  pnt(2,5)=zs(m1,mm2,m)-aa(3)
                  pnt(1,6)=ys(m1+1,mm2,m)-aa(2)
                  pnt(2,6)=zs(m1+1,mm2,m)-aa(3)
                  call point_polygon(6,pnt,keep)
                endif
              endif

              if(keep.eq.1) then
                q1=-(yb*(z(k)-f1(3))-zb*(y(j)-f1(2)))/xn
                q2=-(zt*(y(j)-f1(2))-yt*(z(k)-f1(3)))/xn
                sx=(f2(1)-f1(1))*q1+(f4(1)-f1(1))*q2+f1(1)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sx
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f2(n)-f1(n))*q1+(f4(n)-f1(n))*q2+f1(n)
                enddo

                i=int((sx-x0)/hdx)
                if(sx.eq.x(i).and.sign(1.0,xn).gt.0.0) i=i-1
                if(mod(i,2).eq.0) then
                  lc=i/2+1
                  le=lc-1
                else
                  le=(i+1)/2
                  lc=le
                endif
                if(mod(j,2).ne.0) then
                  if(mod(k,2).eq.0) then
                    jekc(1,njekc(m),m)=(j+1)/2
                    jekc(2,njekc(m),m)=k/2+1
                    jekc(3,njekc(m),m)=le
                    jekc(4,njekc(m),m)=lc
                    do n=1,55
                      fjekc(n,njekc(m),m)=f(n)
                    enddo
                    njekc(m)=njekc(m)+1
                  endif
                else
                  if(mod(k,2).ne.0) then
                    jcke(1,njcke(m),m)=j/2+1
                    jcke(2,njcke(m),m)=(k+1)/2
                    jcke(3,njcke(m),m)=le
                    jcke(4,njcke(m),m)=lc
                    do n=1,55
                      fjcke(n,njcke(m),m)=f(n)
                    enddo
                    njcke(m)=njcke(m)+1
                  else
                    jckc(1,njckc(m),m)=j/2+1
                    jckc(2,njckc(m),m)=k/2+1
                    jckc(3,njckc(m),m)=le
                    jckc(4,njckc(m),m)=lc
                    do n=1,55
                      fjckc(n,njckc(m),m)=f(n)
                    enddo
                    njckc(m)=njckc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

c  A=3, B=4, C=2

          aa(1)=f3(1)
          aa(2)=f3(2)
          aa(3)=f3(3)
          bb(1)=f4(1)
          bb(2)=f4(2)
          bb(3)=f4(3)
          cc(1)=f2(1)
          cc(2)=f2(2)
          cc(3)=f2(3)

          l1=int((aa(1)-x0)/hdx)
          l2=int((bb(1)-x0)/hdx)
          l3=int((cc(1)-x0)/hdx)
          call maxmin(l1,l2,l3,iu,id)

          l1=int((aa(2)-y0)/hdy)
          l2=int((bb(2)-y0)/hdy)
          l3=int((cc(2)-y0)/hdy)
          call maxmin(l1,l2,l3,ju,jd)

          l1=int((aa(3)-z0)/hdz)
          l2=int((bb(3)-z0)/hdz)
          l3=int((cc(3)-z0)/hdz)
          call maxmin(l1,l2,l3,ku,kd)

          xt=aa(1)-bb(1)
          yt=aa(2)-bb(2)
          zt=aa(3)-bb(3)
          xb=aa(1)-cc(1)
          yb=aa(2)-cc(2)
          zb=aa(3)-cc(3)
          xn=yt*zb-zt*yb
          yn=zt*xb-xt*zb
          zn=xt*yb-yt*xb
          sn=dsqrt(xn*xn+yn*yn+zn*zn)

c  ij-line

          IF(dabs(zn)/sn.gt.1.0d-12) THEN
          do j=jd,ju
            do i=id,iu
              call line_side(1,2,bb,aa,cc,x(i),y(j),lab)
              call line_side(1,2,cc,aa,bb,x(i),y(j),lac)
              call line_side(1,2,cc,bb,aa,x(i),y(j),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                if(zn0*(sign(1.0,zn)+sign(1.0,zn0)).ne.0.0) keep=1
              endif

              if(keep.eq.1) then
                q1=(yb*(x(i)+f3(1)-f2(1)-f4(1))-
     .              xb*(y(j)+f3(2)-f2(2)-f4(2)))/zn
                q2=(xt*(y(j)+f3(2)-f2(2)-f4(2))-
     .              yt*(x(i)+f3(1)-f2(1)-f4(1)))/zn
                sz=(f3(3)-f4(3))*q1+(f3(3)-f2(3))*q2+f2(3)+f4(3)-f3(3)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sz
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f3(n)-f4(n))*q1+(f3(n)-f2(n))*q2
     .                 +f2(n)+f4(n)-f3(n)
                enddo

                k=int((sz-z0)/hdz)
                if(sz.eq.z(k).and.sign(1.0,zn).gt.0.0) k=k-1
                if(mod(k,2).eq.0) then
                  lc=k/2+1
                  le=lc-1
                else
                  le=(k+1)/2
                  lc=le
                endif
                if(mod(i,2).ne.0) then
                  if(mod(j,2).eq.0) then
                    iejc(1,niejc(m),m)=(i+1)/2
                    iejc(2,niejc(m),m)=j/2+1
                    iejc(3,niejc(m),m)=le
                    iejc(4,niejc(m),m)=lc
                    do n=1,55
                      fiejc(n,niejc(m),m)=f(n)
                    enddo
                    niejc(m)=niejc(m)+1
                  endif
                else
                  if(mod(j,2).ne.0) then
                    icje(1,nicje(m),m)=i/2+1
                    icje(2,nicje(m),m)=(j+1)/2
                    icje(3,nicje(m),m)=le
                    icje(4,nicje(m),m)=lc
                    do n=1,55
                      ficje(n,nicje(m),m)=f(n)
                    enddo
                    nicje(m)=nicje(m)+1
                  else
                    icjc(1,nicjc(m),m)=i/2+1
                    icjc(2,nicjc(m),m)=j/2+1
                    icjc(3,nicjc(m),m)=le
                    icjc(4,nicjc(m),m)=lc
                    do n=1,55
                      ficjc(n,nicjc(m),m)=f(n)
                    enddo
                    nicjc(m)=nicjc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

c  ik-line

          IF(dabs(yn)/sn.gt.1.0d-12) THEN
          do k=kd,ku
            do i=id,iu
              call line_side(1,3,bb,aa,cc,x(i),z(k),lab)
              call line_side(1,3,cc,aa,bb,x(i),z(k),lac)
              call line_side(1,3,cc,bb,aa,x(i),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                if(yn0*(sign(1.0,yn)+sign(1.0,yn0)).ne.0.0) keep=1
              endif

              if(keep.eq.1) then
                q1=-(zb*(x(i)+f3(1)-f2(1)-f4(1))-
     .               xb*(z(k)+f3(3)-f2(3)-f4(3)))/yn
                q2=-(xt*(z(k)+f3(3)-f2(3)-f4(3))-
     .               zt*(x(i)+f3(1)-f2(1)-f4(1)))/yn
                sy=(f3(2)-f4(2))*q1+(f3(2)-f2(2))*q2+f2(2)+f4(2)-f3(2)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sy
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f3(n)-f4(n))*q1+(f3(n)-f2(n))*q2
     .                 +f2(n)+f4(n)-f3(n)
                enddo

                j=int((sy-y0)/hdy)
                if(sy.eq.y(j).and.sign(1.0,yn).gt.0.0) j=j-1
                if(mod(j,2).eq.0) then
                  lc=j/2+1
                  le=lc-1
                else
                  le=(j+1)/2
                  lc=le
                endif
                if(mod(i,2).ne.0) then
                  if(mod(k,2).eq.0) then
                    iekc(1,niekc(m),m)=(i+1)/2
                    iekc(2,niekc(m),m)=k/2+1
                    iekc(3,niekc(m),m)=le
                    iekc(4,niekc(m),m)=lc
                    do n=1,55
                      fiekc(n,niekc(m),m)=f(n)
                    enddo
                    niekc(m)=niekc(m)+1
                  endif
                else
                  if(mod(k,2).ne.0) then
                    icke(1,nicke(m),m)=i/2+1
                    icke(2,nicke(m),m)=(k+1)/2
                    icke(3,nicke(m),m)=le
                    icke(4,nicke(m),m)=lc
                    do n=1,55
                      ficke(n,nicke(m),m)=f(n)
                    enddo
                    nicke(m)=nicke(m)+1
                  else
                    ickc(1,nickc(m),m)=i/2+1
                    ickc(2,nickc(m),m)=k/2+1
                    ickc(3,nickc(m),m)=le
                    ickc(4,nickc(m),m)=lc
                    do n=1,55
                      fickc(n,nickc(m),m)=f(n)
                    enddo
                    nickc(m)=nickc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

c  jk-line

          IF(dabs(xn)/sn.gt.1.0d-12) THEN
          do k=kd,ku
            do j=jd,ju
              call line_side(2,3,bb,aa,cc,y(j),z(k),lab)
              call line_side(2,3,cc,aa,bb,y(j),z(k),lac)
              call line_side(2,3,cc,bb,aa,y(j),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                if(xn0*(sign(1.0,xn)+sign(1.0,xn0)).ne.0.0) keep=1
              endif

              if(keep.eq.1) then
                q1=-(yb*(z(k)+f3(3)-f2(3)-f4(3))-
     .               zb*(y(j)+f3(2)-f2(2)-f4(2)))/xn
                q2=-(zt*(y(j)+f3(2)-f2(2)-f4(2))-
     .               yt*(z(k)+f3(3)-f2(3)-f4(3)))/xn
                sx=(f3(1)-f4(1))*q1+(f3(1)-f2(1))*q2+f2(1)+f4(1)-f3(1)
                if(isphr.eq.0) then
                  f(1)=alfa10(m1)+q1*dalfa10
                else
                  f(1)=alfa11(m1)+q1*dalfa11
                endif
                f(2)=alfa20(m2)+q2*dalfa20
                f(3)=sx
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=xn
                f(11)=yn
                f(12)=zn
                do n=13,55
                  f(n)=(f3(n)-f4(n))*q1+(f3(n)-f2(n))*q2
     .                 +f2(n)+f4(n)-f3(n)
                enddo

                i=int((sx-x0)/hdx)
                if(sx.eq.x(i).and.sign(1.0,xn).gt.0.0) i=i-1
                if(mod(i,2).eq.0) then
                  lc=i/2+1
                  le=lc-1
                else
                  le=(i+1)/2
                  lc=le
                endif
                if(mod(j,2).ne.0) then
                  if(mod(k,2).eq.0) then
                    jekc(1,njekc(m),m)=(j+1)/2
                    jekc(2,njekc(m),m)=k/2+1
                    jekc(3,njekc(m),m)=le
                    jekc(4,njekc(m),m)=lc
                    do n=1,55
                      fjekc(n,njekc(m),m)=f(n)
                    enddo
                    njekc(m)=njekc(m)+1
                  endif
                else
                  if(mod(k,2).ne.0) then
                    jcke(1,njcke(m),m)=j/2+1
                    jcke(2,njcke(m),m)=(k+1)/2
                    jcke(3,njcke(m),m)=le
                    jcke(4,njcke(m),m)=lc
                    do n=1,55
                      fjcke(n,njcke(m),m)=f(n)
                    enddo
                    njcke(m)=njcke(m)+1
                  else
                    jckc(1,njckc(m),m)=j/2+1
                    jckc(2,njckc(m),m)=k/2+1
                    jckc(3,njckc(m),m)=le
                    jckc(4,njckc(m),m)=lc
                    do n=1,55
                      fjckc(n,njckc(m),m)=f(n)
                    enddo
                    njckc(m)=njckc(m)+1
                  endif
                endif

              endif

            enddo
          enddo
          ENDIF

        enddo
      enddo

      if(isphr.eq.1) then
        call poles_cover(m,0,1,1.0d0)
        call poles_cover(m,ns1,ns1-1,-1.0d0)
      endif

      niejc(m)=niejc(m)-1
      nicje(m)=nicje(m)-1
      nicjc(m)=nicjc(m)-1
      niekc(m)=niekc(m)-1
      nicke(m)=nicke(m)-1
      nickc(m)=nickc(m)-1
      njekc(m)=njekc(m)-1
      njcke(m)=njcke(m)-1
      njckc(m)=njckc(m)-1

      imod(1)=mod(niejc(m),2)
      imod(2)=mod(nicje(m),2)
      imod(3)=mod(nicjc(m),2)
      imod(4)=mod(niekc(m),2)
      imod(5)=mod(nicke(m),2)
      imod(6)=mod(nickc(m),2)
      imod(7)=mod(njekc(m),2)
      imod(8)=mod(njcke(m),2)
      imod(9)=mod(njckc(m),2)

      iprod=1
      do im=1,9
        iprod=iprod*imod(im)
      enddo
      if(iprod.eq.0) then
        write(*,*)'  !! warn: odd number of intersection points !!'
      endif

      ENDDO

      return
      end


c-----------------------------------------------------------------------
