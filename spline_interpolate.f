c-----------------------------------------------------------------------
c
      subroutine spline_interpolate(boo,foo,isphr)
      include 'parameter.inc'
      include 'surface.inc'
      integer n1,n2,nd,nh,isphr
      parameter(n1=ks1-1,n2=ks2-1,nd=2*n1+2,nh=ns2/2-1)
      real*8 foo(0:ns1,0:ns2)
      real*8 alfa0(0:n1),alfa1(0:nd),alfa2(0:n2),alfa3(0:2*ns1+2)
      real*8 boo(0:n1,0:n2),coo(0:n1,0:n2),goo(0:n1,0:ns2)
      real*8 f(0:n1),f1(0:n1),cf(0:n1),cf1(0:n1),ff(0:n1),dff(0:n1)
      real*8 doo(0:nd,0:nh),eoo(0:nd,0:nh)
      real*8 b(0:nh),b1(0:nh),cb(0:nh),cb1(0:nh),bb(0:nh),dbb(0:nh)
      real*8 hoo(0:n1,0:ns2)
      real*8 g(0:ns2),g1(0:ns2),cg(0:ns2),cg1(0:ns2)
      real*8 gg(0:ns2),dgg(0:ns2)

      do k2=0,n2
        m2=k2*kip2
        alfa2(k2)=alfa20(m2)
      enddo
      do k1=0,n1
        m1=k1*kip1
        alfa0(k1)=alfa10(m1)
        alfa1(k1)=alfa11(m1)
        alfa1(n1+1+k1)=2.0d0*sl1-alfa11(ns1-m1)
      enddo
      alfa1(nd)=alfa1(0)+2.0d0*sl1

      call cubic_spline2(alfa2,boo,coo,n2,n2,n1+1)

      do m2=0,ns2-1
        k2=int(m2/kip2)
        do k1=0,n1
          f(k1)=boo(k1,k2)
          f1(k1)=boo(k1,k2+1)
          cf(k1)=coo(k1,k2)
          cf1(k1)=coo(k1,k2+1)
        enddo
        call fdf(alfa2(k2),alfa2(k2+1),f,f1,cf,cf1,
     .           alfa20(m2),ff,dff,n1+1)
        do k1=0,n1
          goo(k1,m2)=ff(k1)
        enddo
      enddo
      do k1=0,n1
        goo(k1,ns2)=goo(k1,0)
      enddo

      if(isphr.eq.1) then

        do m1=0,ns1
          alfa3(m1)=alfa11(m1)
          alfa3(ns1+1+m1)=2.0d0*sl1-alfa11(ns1-m1)
        enddo
        alfa3(2*ns1+2)=alfa3(0)+2.0d0*sl1

        do m2=0,nh
          do k1=0,n1
            doo(k1,m2)=goo(k1,m2)
          enddo
          do k1=n1+1,nd-1
            doo(k1,m2)=goo(nd-1-k1,m2+nh+1)
          enddo
          doo(nd,m2)=doo(0,m2)
        enddo

        call cubic_spline1(alfa1,doo,eoo,nd,nd,nh+1)

        do m1=0,2*ns1+1
          if(m1.le.ns1) then
            k1=int(m1/kip1)
          else
            k1=n1+int((m1-ns1-1)/kip1)+1
          endif
          do m2=0,nh
            b(m2)=doo(k1,m2)
            b1(m2)=doo(k1+1,m2)
            cb(m2)=eoo(k1,m2)
            cb1(m2)=eoo(k1+1,m2)
          enddo
          call fdf(alfa1(k1),alfa1(k1+1),b,b1,cb,cb1,
     .             alfa3(m1),bb,dbb,nh+1)
          if(m1.le.ns1) then
            do m2=0,nh
              foo(m1,m2)=bb(m2)
            enddo
          else
            do m2=0,nh
              foo(2*ns1+1-m1,m2+nh+1)=bb(m2)
            enddo
          endif
        enddo

      endif


      if(isphr.eq.0) then

        call cubic_spline1(alfa0,goo,hoo,n1,n1,ns2+1)

        do m1=0,ns1-1
          k1=int(m1/kip1)
          do m2=0,ns2
            g(m2)=goo(k1,m2)
            g1(m2)=goo(k1+1,m2)
            cg(m2)=hoo(k1,m2)
            cg1(m2)=hoo(k1+1,m2)
          enddo
          call fdf(alfa0(k1),alfa0(k1+1),g,g1,cg,cg1,
     .             alfa10(m1),gg,dgg,ns2+1)
          do m2=0,ns2
            foo(m1,m2)=gg(m2)
          enddo
        enddo
        do m2=0,ns2
          foo(ns1,m2)=foo(0,m2)
        enddo

      endif

      do m1=0,ns1
        foo(m1,ns2)=foo(m1,0)
      enddo

      return
      end


c-----------------------------------------------------------------------
