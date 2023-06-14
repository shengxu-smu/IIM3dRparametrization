c-----------------------------------------------------------------------
c
      subroutine c2_solver(a,b,c,x,y,z,r,v,icheck)
      implicit none
      integer icheck
      real*8 a,b,c,x,y,z
      real*8 r(7),v(6)
      real*8 s2,s3,el,em,en,d2,d3,w,e,f,g,h,v5,v6

      if(a.eq.0.0d0) then
        write(*,*)' !! wrong arrangement in matrix !!'
        stop
      endif

      s2=a*a+c*c
      s3=a*a+b*b
      d2=a*x+c*z
      d3=a*x+b*y
      w=c*y+b*z
      el=b*z-c*y
      em=c*x-a*z
      en=a*y-b*x
      if((el*el+em*em+en*en).eq.0.0d0) then
        write(*,*)' !! parallel vectors !!'
        stop
      endif
      e=-w*s3+2.0d0*b*c*d3
      f=-d2*s3+s2*d3
      g=-em*s3-2.0d0*b*c*en
      h=-s2*en

      r(2)=r(2)-a*r(1)
      r(3)=r(3)-x*r(1)

      r(2)=a*r(2)-b*r(4)
      r(3)=a*r(3)-y*r(4)
      r(5)=a*r(5)-x*r(4)

      r(2)=r(2)-c*r(6)
      r(3)=r(3)-z*r(6)
      r(7)=a*r(7)-x*r(6)

      r(3)=s3*r(3)-d3*r(2)
      r(5)=s3*r(5)+en*r(2)

      if(abs(en).ge.abs(em)) then
        r(3)=en*r(3)-e*r(7)
        r(5)=en*r(5)-g*r(7)
        v(6)=r(5)/(en*h+em*g)
        v(5)=(r(7)+em*v(6))/en
        if(icheck.eq.1) then
          if((en*f+em*e).ne.0.0d0) then
            write(*,*)v(6),r(5)/(en*h+em*g)
          else
            write(*,*)en*f+em*e,r(3)
          endif
        endif
      else
        r(3)=em*r(3)+f*r(7)
        r(5)=em*r(5)+h*r(7)
        v(5)=r(5)/(en*h+em*g)
        v(6)=-(r(7)-en*v(5))/em
        if(icheck.eq.1) then
          if((en*f+em*e).ne.0.0d0) then
            write(*,*)v(5),r(5)/(en*h+em*g)
          else
            write(*,*)en*f+em*e,r(3)
          endif
        endif
      endif

      v(4)=-(r(2)+2.0d0*b*c*v(5)+s2*v(6))/s3
      v(3)=(r(6)-b*v(5)-c*v(6))/a
      v(2)=(r(4)-b*v(4)-c*v(5))/a
      v(1)=r(1)-v(4)-v(6)

      return
      end


c-----------------------------------------------------------------------

