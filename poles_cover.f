c-----------------------------------------------------------------------
c
      subroutine poles_cover(m,m1,mm1,flip)
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer l1,l2,l3,iu,id,ju,jd,ku,kd,lc,le
      integer mm1,mm2,lab,lac,lbc,keep,nij,nik,njk,nlimit
      real*8 flip,xt,yt,zt,xb,yb,zb,xn,yn,zn
      real*8 sx,sy,sz,tx,ty,tz,bx,by,bz,an,bn,cn,wa,wb,wc,wsum
      real*8 aa(3),bb(3),cc(3),f1(55),f2(55),f3(55),f(55)
      real*8 pnt0(2,ns2+1),pnt1(2,4),pnt(2,5)

      nij=0
      nik=0
      njk=0
      nlimit=2

      do m2=1,ns2-2

          aa(1)=xs(m1,0,m)
          bb(1)=xs(m1,m2,m)
          cc(1)=xs(m1,m2+1,m)

          aa(2)=ys(m1,0,m)
          bb(2)=ys(m1,m2,m)
          cc(2)=ys(m1,m2+1,m)

          aa(3)=zs(m1,0,m)
          bb(3)=zs(m1,m2,m)
          cc(3)=zs(m1,m2+1,m)

          do n=1,3
            f1(n+12)=ujc1(m1,0,n,m)
            f2(n+12)=ujc1(m1,m2,n,m)
            f3(n+12)=ujc1(m1,m2+1,n,m)

            f1(n+15)=vjc1(m1,0,n,m)
            f2(n+15)=vjc1(m1,m2,n,m)
            f3(n+15)=vjc1(m1,m2+1,n,m)

            f1(n+18)=wjc1(m1,0,n,m)
            f2(n+18)=wjc1(m1,m2,n,m)
            f3(n+18)=wjc1(m1,m2+1,n,m)

            f1(n+21)=pjc1(m1,0,n,m)
            f2(n+21)=pjc1(m1,m2,n,m)
            f3(n+21)=pjc1(m1,m2+1,n,m)
          enddo

          do n=1,6
            f1(n+24)=ujc2(m1,0,n,m)
            f2(n+24)=ujc2(m1,m2,n,m)
            f3(n+24)=ujc2(m1,m2+1,n,m)

            f1(n+30)=vjc2(m1,0,n,m)
            f2(n+30)=vjc2(m1,m2,n,m)
            f3(n+30)=vjc2(m1,m2+1,n,m)

            f1(n+36)=wjc2(m1,0,n,m)
            f2(n+36)=wjc2(m1,m2,n,m)
            f3(n+36)=wjc2(m1,m2+1,n,m)

            f1(n+42)=pjc2(m1,0,n,m)
            f2(n+42)=pjc2(m1,m2,n,m)
            f3(n+42)=pjc2(m1,m2+1,n,m)

            f1(n+49)=pjc(m1,0,n,m)
            f2(n+49)=pjc(m1,m2,n,m)
            f3(n+49)=pjc(m1,m2+1,n,m)
          enddo

          f1(49)=pjc0(m1,0,m)
          f2(49)=pjc0(m1,m2,m)
          f3(49)=pjc0(m1,m2+1,m)

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

          xt=cc(1)-aa(1)
          yt=cc(2)-aa(2)
          zt=cc(3)-aa(3)
          xb=cc(1)-bb(1)
          yb=cc(2)-bb(2)
          zb=cc(3)-bb(3)
          xn=yt*zb-zt*yb
          yn=zt*xb-xt*zb
          zn=xt*yb-yt*xb

c  ij-line

          IF(zn.ne.0.0d0) THEN
          do j=jd,ju
            do i=id,iu
              if(m2.eq.1) then
                call line_side(1,2,aa,bb,cc,x(i),y(j),lab)
              else
                call line_side(1,2,bb,aa,cc,x(i),y(j),lab)
              endif
              call line_side(1,2,cc,aa,bb,x(i),y(j),lac)
              call line_side(1,2,bb,cc,aa,x(i),y(j),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  mm2=m2-1
                  bx=bb(1)-aa(1)
                  by=bb(2)-aa(2)
                  tx=xs(mm1,mm2,m)-aa(1)
                  ty=ys(mm1,mm2,m)-aa(2)
                  cn=tx*by-ty*bx
                  if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.eq.ns2-2) then
                  mm2=ns2-1
                  bx=xs(mm1,mm2,m)-cc(1)
                  by=ys(mm1,mm2,m)-cc(2)
                  cn=xt*by-yt*bx
                  if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
                else
                  mm2=m2+2
                  bx=xs(m1,mm2,m)-cc(1)
                  by=ys(m1,mm2,m)-cc(2)
                  cn=xt*by-yt*bx
                  if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                tx=xs(mm1,m2,m)-bb(1)
                ty=ys(mm1,m2,m)-bb(2)
                cn=tx*yb-ty*xb
                if(cn*(sign(1.0,zn)+sign(1.0,cn)).ne.0.0) keep=1
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  do mm2=1,ns2-1
                    pnt0(1,mm2)=xs(m1,mm2,m)-aa(1)
                    pnt0(2,mm2)=ys(m1,mm2,m)-aa(2)
                  enddo
                  pnt0(1,ns2)=xs(mm1,ns2-1,m)-aa(1)
                  pnt0(2,ns2)=ys(mm1,ns2-1,m)-aa(2)
                  pnt0(1,ns2+1)=xs(mm1,ns2,m)-aa(1)
                  pnt0(2,ns2+1)=ys(mm1,ns2,m)-aa(2)
                  call point_polygon(ns2+1,pnt0,keep)
                endif
              endif
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.0) then
                if(m2.ne.1) then
                  keep=0
                else
                  pnt1(1,1)=xs(m1,m2+1,m)-bb(1)
                  pnt1(2,1)=ys(m1,m2+1,m)-bb(2)
                  pnt1(1,2)=xs(m1,m2-1,m)-bb(1)
                  pnt1(2,2)=ys(m1,m2-1,m)-bb(2)
                  pnt1(1,3)=xs(mm1,m2-1,m)-bb(1)
                  pnt1(2,3)=ys(mm1,m2-1,m)-bb(2)
                  pnt1(1,4)=xs(mm1,m2,m)-bb(1)
                  pnt1(2,4)=ys(mm1,m2,m)-bb(2)
                  call point_polygon(4,pnt1,keep)
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.0) then
                if(m2.ne.ns2-2) then
                  pnt(1,1)=xs(m1,0,m)-cc(1)
                  pnt(2,1)=ys(m1,0,m)-cc(2)
                  pnt(1,2)=xs(m1,m2,m)-cc(1)
                  pnt(2,2)=ys(m1,m2,m)-cc(2)
                  pnt(1,3)=xs(mm1,m2,m)-cc(1)
                  pnt(2,3)=ys(mm1,m2,m)-cc(2)
                  pnt(1,4)=xs(mm1,m2+1,m)-cc(1)
                  pnt(2,4)=ys(mm1,m2+1,m)-cc(2)
                  pnt(1,5)=xs(m1,m2+2,m)-cc(1)
                  pnt(2,5)=ys(m1,m2+2,m)-cc(2)
                  call point_polygon(5,pnt,keep)
                else
                  pnt1(1,1)=xs(m1,0,m)-cc(1)
                  pnt1(2,1)=ys(m1,0,m)-cc(2)
                  pnt1(1,2)=xs(m1,m2,m)-cc(1)
                  pnt1(2,2)=ys(m1,m2,m)-cc(2)
                  pnt1(1,3)=xs(mm1,m2,m)-cc(1)
                  pnt1(2,3)=ys(mm1,m2,m)-cc(2)
                  pnt1(1,4)=xs(mm1,m2+1,m)-cc(1)
                  pnt1(2,4)=ys(mm1,m2+1,m)-cc(2)
                  call point_polygon(4,pnt1,keep)
                endif
              endif

              if(keep.eq.1) then
                nij=nij+1
                tx=bb(1)-x(i)
                ty=bb(2)-y(j)
                bx=cc(1)-x(i)
                by=cc(2)-y(j)
                wa=dabs(tx*by-ty*bx)
                tx=cc(1)-x(i)
                ty=cc(2)-y(j)
                bx=aa(1)-x(i)
                by=aa(2)-y(j)
                wb=dabs(tx*by-ty*bx)
                tx=aa(1)-x(i)
                ty=aa(2)-y(j)
                bx=bb(1)-x(i)
                by=bb(2)-y(j)
                wc=dabs(tx*by-ty*bx)
                wsum=wa+wb+wc
                wa=wa/wsum
                wb=wb/wsum
                wc=wc/wsum
                sz=wa*aa(3)+wb*bb(3)+wc*cc(3)
                f(1)=0.5d0*sl1*(1.0d0-flip)
                f(2)=wa*alfa20(0)+wb*alfa20(m2)+wc*alfa20(m2+1)
                f(3)=sz
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=flip*xn
                f(11)=flip*yn
                f(12)=flip*zn
                do n=13,55
                  f(n)=wa*f1(n)+wb*f2(n)+wc*f3(n)
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
c                   ijdex(i/2+1,j/2+1,lc)=nicjc(m)
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

          IF(yn.ne.0.0d0) THEN
          do k=kd,ku
            do i=id,iu
              if(m2.eq.1) then
                call line_side(1,3,aa,bb,cc,x(i),z(k),lab)
              else
                call line_side(1,3,bb,aa,cc,x(i),z(k),lab)
              endif
              call line_side(1,3,cc,aa,bb,x(i),z(k),lac)
              call line_side(1,3,bb,cc,aa,x(i),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  mm2=m2-1
                  bx=bb(1)-aa(1)
                  bz=bb(3)-aa(3)
                  tx=xs(mm1,mm2,m)-aa(1)
                  tz=zs(mm1,mm2,m)-aa(3)
                  bn=tz*bx-tx*bz
                  if(bn*(sign(1.0,yn)+sign(1.0,bn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.eq.ns2-2) then
                  mm2=ns2-1
                  bx=xs(mm1,mm2,m)-cc(1)
                  bz=zs(mm1,mm2,m)-cc(3)
                  bn=zt*bx-xt*bz
                  if(bn*(sign(1.0,yn)+sign(1.0,bn)).ne.0.0) keep=1
                else
                  mm2=m2+2
                  bx=xs(m1,mm2,m)-cc(1)
                  bz=zs(m1,mm2,m)-cc(3)
                  bn=zt*bx-xt*bz
                  if(bn*(sign(1.0,bn)+sign(1.0,bn)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                tx=xs(mm1,m2,m)-bb(1)
                tz=zs(mm1,m2,m)-bb(3)
                bn=tz*xb-tx*zb
                if(bn*(sign(1.0,yn)+sign(1.0,bn)).ne.0.0) keep=1
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  do mm2=1,ns2-1
                    pnt0(1,mm2)=xs(m1,mm2,m)-aa(1)
                    pnt0(2,mm2)=zs(m1,mm2,m)-aa(3)
                  enddo
                  pnt0(1,ns2)=xs(mm1,ns2-1,m)-aa(1)
                  pnt0(2,ns2)=zs(mm1,ns2-1,m)-aa(3)
                  pnt0(1,ns2+1)=xs(mm1,ns2,m)-aa(1)
                  pnt0(2,ns2+1)=zs(mm1,ns2,m)-aa(3)
                  call point_polygon(ns2+1,pnt0,keep)
                endif
              endif
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.0) then
                if(m2.ne.1) then
                  keep=0
                else
                  pnt1(1,1)=xs(m1,m2+1,m)-bb(1)
                  pnt1(2,1)=zs(m1,m2+1,m)-bb(3)
                  pnt1(1,2)=xs(m1,m2-1,m)-bb(1)
                  pnt1(2,2)=zs(m1,m2-1,m)-bb(3)
                  pnt1(1,3)=xs(mm1,m2-1,m)-bb(1)
                  pnt1(2,3)=zs(mm1,m2-1,m)-bb(3)
                  pnt1(1,4)=xs(mm1,m2,m)-bb(1)
                  pnt1(2,4)=zs(mm1,m2,m)-bb(3)
                  call point_polygon(4,pnt1,keep)
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.0) then
                if(m2.ne.ns2-2) then
                  pnt(1,1)=xs(m1,0,m)-cc(1)
                  pnt(2,1)=zs(m1,0,m)-cc(3)
                  pnt(1,2)=xs(m1,m2,m)-cc(1)
                  pnt(2,2)=zs(m1,m2,m)-cc(3)
                  pnt(1,3)=xs(mm1,m2,m)-cc(1)
                  pnt(2,3)=zs(mm1,m2,m)-cc(3)
                  pnt(1,4)=xs(mm1,m2+1,m)-cc(1)
                  pnt(2,4)=zs(mm1,m2+1,m)-cc(3)
                  pnt(1,5)=xs(m1,m2+2,m)-cc(1)
                  pnt(2,5)=zs(m1,m2+2,m)-cc(3)
                  call point_polygon(5,pnt,keep)
                else
                  pnt1(1,1)=xs(m1,0,m)-cc(1)
                  pnt1(2,1)=zs(m1,0,m)-cc(3)
                  pnt1(1,2)=xs(m1,m2,m)-cc(1)
                  pnt1(2,2)=zs(m1,m2,m)-cc(3)
                  pnt1(1,3)=xs(mm1,m2,m)-cc(1)
                  pnt1(2,3)=zs(mm1,m2,m)-cc(3)
                  pnt1(1,4)=xs(mm1,m2+1,m)-cc(1)
                  pnt1(2,4)=zs(mm1,m2+1,m)-cc(3)
                  call point_polygon(4,pnt1,keep)
                endif
              endif

              if(keep.eq.1) then
                nik=nik+1
                tx=bb(1)-x(i)
                tz=bb(3)-z(k)
                bx=cc(1)-x(i)
                bz=cc(3)-z(k)
                wa=dabs(tx*bz-tz*bx)
                tx=cc(1)-x(i)
                tz=cc(3)-z(k)
                bx=aa(1)-x(i)
                bz=aa(3)-z(k)
                wb=dabs(tx*bz-tz*bx)
                tx=aa(1)-x(i)
                tz=aa(3)-z(j)
                bx=bb(1)-x(i)
                bz=bb(3)-z(j)
                wc=dabs(tx*bz-tz*bx)
                wsum=wa+wb+wc
                wa=wa/wsum
                wb=wb/wsum
                wc=wc/wsum
                sy=wa*aa(2)+wb*bb(2)+wc*cc(2)
                f(1)=0.5d0*sl1*(1.0d0-flip)
                f(2)=wa*alfa20(0)+wb*alfa20(m2)+wc*alfa20(m2+1)
                f(3)=sy
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=flip*xn
                f(11)=flip*yn
                f(12)=flip*zn
                do n=13,55
                  f(n)=wa*f1(n)+wb*f2(n)+wc*f3(n)
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
c                   ikdex(i/2+1,lc,k/2+1)=nickc(m)
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

          IF(xn.ne.0.0d0) THEN
          do k=kd,ku
            do j=jd,ju
              if(m2.eq.1) then
                call line_side(2,3,aa,bb,cc,y(j),z(k),lab)
              else
                call line_side(2,3,bb,aa,cc,y(j),z(k),lab)
              endif
              call line_side(2,3,cc,aa,bb,y(j),z(k),lac)
              call line_side(2,3,bb,cc,aa,y(j),z(k),lbc)
              keep=0
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.1) keep=1
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  mm2=m2-1
                  by=bb(2)-aa(2)
                  bz=bb(3)-aa(3)
                  ty=ys(mm1,mm2,m)-aa(2)
                  tz=zs(mm1,mm2,m)-aa(3)
                  an=ty*bz-tz*by
                  if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.eq.ns2-2) then
                  mm2=ns2-1
                  by=ys(mm1,mm2,m)-cc(2)
                  bz=zs(mm1,mm2,m)-cc(3)
                  an=yt*bz-zt*by
                  if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
                else
                  mm2=m2+2
                  by=ys(m1,mm2,m)-cc(2)
                  bz=zs(m1,mm2,m)-cc(3)
                  an=yt*bz-zt*by
                  if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
                endif
              endif
              if(lab.eq.1.and.lac.eq.1.and.lbc.eq.0) then
                ty=ys(mm1,m2,m)-bb(2)
                tz=zs(mm1,m2,m)-bb(3)
                an=ty*zb-tz*yb
                if(an*(sign(1.0,xn)+sign(1.0,an)).ne.0.0) keep=1
              endif
              if(lab.eq.0.and.lac.eq.0.and.lbc.eq.1) then
                if(m2.ne.1) then
                  keep=0
                else
                  do mm2=1,ns2-1
                    pnt0(1,mm2)=ys(m1,mm2,m)-aa(2)
                    pnt0(2,mm2)=zs(m1,mm2,m)-aa(3)
                  enddo
                  pnt0(1,ns2)=ys(mm1,ns2-1,m)-aa(2)
                  pnt0(2,ns2)=zs(mm1,ns2-1,m)-aa(3)
                  pnt0(1,ns2+1)=ys(mm1,ns2,m)-aa(2)
                  pnt0(2,ns2+1)=zs(mm1,ns2,m)-aa(3)
                  call point_polygon(ns2+1,pnt0,keep)
                endif
              endif
              if(lab.eq.0.and.lac.eq.1.and.lbc.eq.0) then
                if(m2.ne.1) then
                  keep=0
                else
                  pnt1(1,1)=ys(m1,m2+1,m)-bb(2)
                  pnt1(2,1)=zs(m1,m2+1,m)-bb(3)
                  pnt1(1,2)=ys(m1,m2-1,m)-bb(2)
                  pnt1(2,2)=zs(m1,m2-1,m)-bb(3)
                  pnt1(1,3)=ys(mm1,m2-1,m)-bb(2)
                  pnt1(2,3)=zs(mm1,m2-1,m)-bb(3)
                  pnt1(1,4)=ys(mm1,m2,m)-bb(2)
                  pnt1(2,4)=zs(mm1,m2,m)-bb(3)
                  call point_polygon(4,pnt1,keep)
                endif
              endif
              if(lab.eq.1.and.lac.eq.0.and.lbc.eq.0) then
                if(m2.ne.ns2-2) then
                  pnt(1,1)=ys(m1,0,m)-cc(2)
                  pnt(2,1)=zs(m1,0,m)-cc(3)
                  pnt(1,2)=ys(m1,m2,m)-cc(2)
                  pnt(2,2)=zs(m1,m2,m)-cc(3)
                  pnt(1,3)=ys(mm1,m2,m)-cc(2)
                  pnt(2,3)=zs(mm1,m2,m)-cc(3)
                  pnt(1,4)=ys(mm1,m2+1,m)-cc(2)
                  pnt(2,4)=zs(mm1,m2+1,m)-cc(3)
                  pnt(1,5)=ys(m1,m2+2,m)-cc(2)
                  pnt(2,5)=zs(m1,m2+2,m)-cc(3)
                  call point_polygon(5,pnt,keep)
                else
                  pnt1(1,1)=ys(m1,0,m)-cc(2)
                  pnt1(2,1)=zs(m1,0,m)-cc(3)
                  pnt1(1,2)=ys(m1,m2,m)-cc(2)
                  pnt1(2,2)=zs(m1,m2,m)-cc(3)
                  pnt1(1,3)=ys(mm1,m2,m)-cc(2)
                  pnt1(2,3)=zs(mm1,m2,m)-cc(3)
                  pnt1(1,4)=ys(mm1,m2+1,m)-cc(2)
                  pnt1(2,4)=zs(mm1,m2+1,m)-cc(3)
                  call point_polygon(4,pnt1,keep)
                endif
              endif

              if(keep.eq.1) then
                njk=njk+1
                ty=bb(2)-y(j)
                tz=bb(3)-z(k)
                by=cc(2)-y(j)
                bz=cc(3)-z(k)
                wa=dabs(ty*bz-tz*by)
                ty=cc(2)-y(j)
                tz=cc(3)-z(k)
                by=aa(2)-y(j)
                bz=aa(3)-z(k)
                wb=dabs(ty*bz-tz*by)
                ty=aa(2)-y(j)
                tz=aa(3)-z(k)
                by=bb(2)-y(j)
                bz=bb(3)-z(k)
                wc=dabs(ty*bz-tz*by)
                wsum=wa+wb+wc
                wa=wa/wsum
                wb=wb/wsum
                wc=wc/wsum
                sx=wa*aa(1)+wb*bb(1)+wc*cc(1)
                f(1)=0.5d0*sl1*(1.0d0-flip)
                f(2)=wa*alfa20(0)+wb*alfa20(m2)+wc*alfa20(m2+1)
                f(3)=sx
                do n=4,9
                  f(n)=0.0d0
                enddo
                f(10)=flip*xn
                f(11)=flip*yn
                f(12)=flip*zn
                do n=13,55
                  f(n)=wa*f1(n)+wb*f2(n)+wc*f3(n)
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
c                   jkdex(lc,j/2+1,k/2+1)=njckc(m)
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

      if(nij.gt.nlimit.or.nik.gt.nlimit.or.njk.gt.nlimit) then
c       write(*,*)'  !! warn: more than 2 point in holes at poles !!'
      endif

      return
      end


c-----------------------------------------------------------------------
