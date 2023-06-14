c-----------------------------------------------------------------------
c
      subroutine index_initialize
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer isphr,inside,l1,l2,l3,la,lb,lc
      integer iu,id,ju,jd,ku,kd,ie,ic,je,jc,ke,kc
      integer n1,n2,n3,n4,n5,n6,n7,n8
      integer io(0:2*nx-2,0:2*ny-2,0:2*nz-2)
      real*8 f1(3),f2(3),f3(3),f4(3),dd(3)

      do k=0,2*nz-2
        do j=0,2*ny-2
          do i=0,2*nx-2
            io(i,j,k)=0
          enddo
        enddo
      enddo

      DO m=1,ms

      isphr=i01(m)

      dd(1)=0.0d0
      dd(2)=0.0d0
      dd(3)=0.0d0
      l1=int((dd(1)-x0)/hdx)
      l2=int((dd(2)-y0)/hdy)
      l3=int((dd(3)-z0)/hdz)

      do m2=0,ns2-1
        do m1=0,ns1-1

          f1(1)=xss(m1,m2,m)
          f2(1)=xss(m1+1,m2,m)
          f3(1)=xss(m1+1,m2+1,m)
          f4(1)=xss(m1,m2+1,m)

          f1(2)=yss(m1,m2,m)
          f2(2)=yss(m1+1,m2,m)
          f3(2)=yss(m1+1,m2+1,m)
          f4(2)=yss(m1,m2+1,m)

          f1(3)=zss(m1,m2,m)
          f2(3)=zss(m1+1,m2,m)
          f3(3)=zss(m1+1,m2+1,m)
          f4(3)=zss(m1,m2+1,m)

c  A=1, B=2, C=4:

          la=int((f1(1)-x0)/hdx)
          lb=int((f2(1)-x0)/hdx)
          lc=int((f4(1)-x0)/hdx)
          call maxmin(la,lb,lc,iu,id)
          if(l1.gt.iu) iu=l1
          if(l1.lt.id) id=l1

          la=int((f1(2)-y0)/hdy)
          lb=int((f2(2)-y0)/hdy)
          lc=int((f4(2)-y0)/hdy)
          call maxmin(la,lb,lc,ju,jd)
          if(l2.gt.ju) ju=l2
          if(l2.lt.jd) jd=l2
          
          la=int((f1(3)-z0)/hdz)
          lb=int((f2(3)-z0)/hdz)
          lc=int((f4(3)-z0)/hdz)
          call maxmin(la,lb,lc,ku,kd)
          if(l3.gt.ku) ku=l3
          if(l3.lt.kd) kd=l3

          IF(isphr.eq.1) then
          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call inside_test(f1,f2,f4,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m
              enddo
            enddo
          enddo
          ENDIF

          IF(isphr.eq.0) then
          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call in124_test(f1,f2,f4,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m-io(i,j,k)
              enddo
            enddo
          enddo
          ENDIF


c  A=3, B=4, C=2

          la=int((f3(1)-x0)/hdx)
          lb=int((f4(1)-x0)/hdx)
          lc=int((f2(1)-x0)/hdx)
          call maxmin(la,lb,lc,iu,id)
          if(l1.gt.iu) iu=l1
          if(l1.lt.id) id=l1

          la=int((f3(2)-y0)/hdy)
          lb=int((f4(2)-y0)/hdy)
          lc=int((f2(2)-y0)/hdy)
          call maxmin(la,lb,lc,ju,jd)
          if(l2.gt.ju) ju=l2
          if(l2.lt.jd) jd=l2

          la=int((f3(3)-z0)/hdz)
          lb=int((f4(3)-z0)/hdz)
          lc=int((f2(3)-z0)/hdz)
          call maxmin(la,lb,lc,ku,kd)
          if(l3.gt.ku) ku=l3
          if(l3.lt.kd) kd=l3

          IF(isphr.eq.1) then
          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call inside_test(f3,f4,f2,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m
              enddo
            enddo
          enddo
          ENDIF

          IF(isphr.eq.0) then
          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call in324_test(f3,f2,f4,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m-io(i,j,k)
              enddo
            enddo
          enddo
          ENDIF

        enddo
      enddo

      IF(isphr.eq.1) THEN

        do m2=1,ns2-2

c  North pole

          f1(1)=xss(0,0,m)
          f2(1)=xss(0,m2,m)
          f4(1)=xss(0,m2+1,m)

          f1(2)=yss(0,0,m)
          f2(2)=yss(0,m2,m)
          f4(2)=yss(0,m2+1,m)

          f1(3)=zss(0,0,m)
          f2(3)=zss(0,m2,m)
          f4(3)=zss(0,m2+1,m)

          la=int((f1(1)-x0)/hdx)
          lb=int((f2(1)-x0)/hdx)
          lc=int((f4(1)-x0)/hdx)
          call maxmin(la,lb,lc,iu,id)
          if(l1.gt.iu) iu=l1
          if(l1.lt.id) id=l1

          la=int((f1(2)-y0)/hdy)
          lb=int((f2(2)-y0)/hdy)
          lc=int((f4(2)-y0)/hdy)
          call maxmin(la,lb,lc,ju,jd)
          if(l2.gt.ju) ju=l2
          if(l2.lt.jd) jd=l2

          la=int((f1(3)-z0)/hdz)
          lb=int((f2(3)-z0)/hdz)
          lc=int((f4(3)-z0)/hdz)
          call maxmin(la,lb,lc,ku,kd)
          if(l3.gt.ku) ku=l3
          if(l3.lt.kd) kd=l3

          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call inside_test(f1,f2,f4,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m
              enddo
            enddo
          enddo

c  South pole

          f3(1)=xss(ns1,0,m)
          f4(1)=xss(ns1,m2,m)
          f2(1)=xss(ns1,m2+1,m)

          f3(2)=yss(ns1,0,m)
          f4(2)=yss(ns1,m2,m)
          f2(2)=yss(ns1,m2+1,m)

          f3(3)=zss(ns1,0,m)
          f4(3)=zss(ns1,m2,m)
          f2(3)=zss(ns1,m2+1,m)

          la=int((f3(1)-x0)/hdx)
          lb=int((f4(1)-x0)/hdx)
          lc=int((f2(1)-x0)/hdx)
          call maxmin(la,lb,lc,iu,id)
          if(l1.gt.iu) iu=l1
          if(l1.lt.id) id=l1

          la=int((f3(2)-y0)/hdy)
          lb=int((f4(2)-y0)/hdy)
          lc=int((f2(2)-y0)/hdy)
          call maxmin(la,lb,lc,ju,jd)
          if(l2.gt.ju) ju=l2
          if(l2.lt.jd) jd=l2

          la=int((f3(3)-z0)/hdz)
          lb=int((f4(3)-z0)/hdz)
          lc=int((f2(3)-z0)/hdz)
          call maxmin(la,lb,lc,ku,kd)
          if(l3.gt.ku) ku=l3
          if(l3.lt.kd) kd=l3

          do k=kd,ku
            do j=jd,ju
              do i=id,iu
                call inside_test(f3,f4,f2,dd,x(i),y(j),z(k),inside)
                if(inside.eq.1) io(i,j,k)=m
              enddo
            enddo
          enddo

        enddo

      ENDIF

      ENDDO

      do k=0,2*nz-3
        do j=0,2*ny-3
          do i=0,2*nx-3
            n1=io(i,j,k)
            n2=io(i+1,j,k)
            n3=io(i+1,j+1,k)
            n4=io(i,j+1,k)
            n5=io(i,j,k+1)
            n6=io(i+1,j,k+1)
            n7=io(i+1,j+1,k+1)
            n8=io(i,j+1,k+1)
            ind(i,j,k)=n1
            if(n1.eq.0) then
              if(n2.ne.0) ind(i,j,k)=-n2
              if(n3.ne.0) ind(i,j,k)=-n3
              if(n4.ne.0) ind(i,j,k)=-n4
              if(n5.ne.0) ind(i,j,k)=-n5
              if(n6.ne.0) ind(i,j,k)=-n6
              if(n7.ne.0) ind(i,j,k)=-n7
              if(n8.ne.0) ind(i,j,k)=-n8
            else
              if((n2*n3*n4*n5*n6*n7*n8).eq.0) ind(i,j,k)=-n1
            endif
          enddo
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
