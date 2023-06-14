c-----------------------------------------------------------------------
c
      subroutine surface_derivative
      include 'parameter.inc'
      include 'surface.inc'
      integer isphr
      real*8 xt,yt,zt,xb,yb,zb,xn,yn,zn,gacobi
      real*8 foo(0:ns1,0:ns2),boo(0:ns1,0:ns2)

      DO m=1,ms

      isphr=i01(m)

c  \frac{\partial X_i}{\partial \alfa_1}, i=x,y,z

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=xs(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,1,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,1,m)=foo(m1,m2)
        enddo
      enddo

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=ys(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,1,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,2,m)=foo(m1,m2)
        enddo
      enddo

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=zs(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,1,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,3,m)=foo(m1,m2)
        enddo
      enddo

c  \frac{\partial X_i}{\partial \alfa_2}, i=x,y,z

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=xs(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,4,m)=foo(m1,m2)
        enddo
      enddo

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=ys(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,5,m)=foo(m1,m2)
        enddo
      enddo

      do m2=0,ns2
        do m1=0,ns1
          boo(m1,m2)=zs(m1,m2,m)
        enddo
      enddo
      call spline_derivative(boo,foo,isphr,2,1)
      do m2=0,ns2
        do m1=0,ns1
          fst(m1,m2,6,m)=foo(m1,m2)
        enddo
      enddo

c  n_x, n_y, n_z

      do m2=0,ns2
        do m1=0,ns1
          xt=fst(m1,m2,1,m)
          yt=fst(m1,m2,2,m)
          zt=fst(m1,m2,3,m)
          xb=fst(m1,m2,4,m)
          yb=fst(m1,m2,5,m)
          zb=fst(m1,m2,6,m)
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
        enddo
      enddo

c  \frac{\partial^2 X_i}{\partial \alfa_1\alfa_1}, i=x,y,z

      do n=1,3
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fst(m1,m2,n,m)
          enddo
        enddo     
        call spline_derivative(boo,foo,isphr,1,2)
        do m2=0,ns2
          do m1=0,ns1
            snd(m1,m2,n,m)=foo(m1,m2)
          enddo
        enddo
      enddo

c  \frac{\partial^2 X_i}{\partial \alfa_2\alfa_2}, i=x,y,z

      do n=4,6
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fst(m1,m2,n,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,2,2)
        do m2=0,ns2
          do m1=0,ns1
            snd(m1,m2,n,m)=foo(m1,m2)
          enddo
        enddo
      enddo

c  \frac{\partial^2 X_i}{\partial \alfa_1\alfa_2}, i=x,y,z

      do n=1,3
        do m2=0,ns2
           do m1=0,ns1
            boo(m1,m2)=fst(m1,m2,n,m)
          enddo
        enddo
        call spline_derivative(boo,foo,isphr,2,1)
        do m2=0,ns2
          do m1=0,ns1
            snd(m1,m2,n+6,m)=foo(m1,m2)
          enddo
        enddo
      enddo

c  \frac{N_i}{\partial \alfa_1} and \frac{N_i}{\partial \alfa_2}, i=x,y,z

      do m2=0,ns2
        do m1=0,ns1
          xt=fst(m1,m2,1,m)
          yt=fst(m1,m2,2,m)
          zt=fst(m1,m2,3,m)
          xb=fst(m1,m2,4,m)
          yb=fst(m1,m2,5,m)
          zb=fst(m1,m2,6,m)
          snd(m1,m2,10,m)=yt*snd(m1,m2,9,m)+zb*snd(m1,m2,2,m)-
     .                    zt*snd(m1,m2,8,m)-yb*snd(m1,m2,3,m)
          snd(m1,m2,11,m)=zt*snd(m1,m2,7,m)+xb*snd(m1,m2,3,m)-
     .                    xt*snd(m1,m2,9,m)-zb*snd(m1,m2,1,m)
          snd(m1,m2,12,m)=xt*snd(m1,m2,8,m)+yb*snd(m1,m2,1,m)-
     .                    yt*snd(m1,m2,7,m)-xb*snd(m1,m2,2,m)
          snd(m1,m2,13,m)=yt*snd(m1,m2,6,m)+zb*snd(m1,m2,8,m)-
     .                    zt*snd(m1,m2,5,m)-yb*snd(m1,m2,9,m)
          snd(m1,m2,14,m)=zt*snd(m1,m2,4,m)+xb*snd(m1,m2,9,m)-
     .                    xt*snd(m1,m2,6,m)-zb*snd(m1,m2,7,m)
          snd(m1,m2,15,m)=xt*snd(m1,m2,5,m)+yb*snd(m1,m2,7,m)-
     .                    yt*snd(m1,m2,4,m)-xb*snd(m1,m2,8,m)
        enddo
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------
