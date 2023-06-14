c-----------------------------------------------------------------------
c
      subroutine qcorrection_difference
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc,ic1,jc1,kc1
      real*8 sx,sy,sz,signx,signy,signz

      DO m=1,ms

      do n=0,nicjc(m)
        signz=sign(1.0,ficjc(12,n,m))
        sz=ficjc(3,n,m)
        ke=icjc(3,n,m)
        kc=icjc(4,n,m)
        kc1=kc+1
        if(ke.eq.kc) then
          pdz(n,m)=-dz1*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc1)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc1)-sz)**2.0d0)
        else
          pdz(n,m)=-dz1*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc)-sz)**2.0d0)
        endif
        pdzz(1,n,m)=-dz2*(ficjc(49,n,m)*signz+
     .                    ficjc(24,n,m)*(zc(kc1)-sz)+
     .              0.5d0*ficjc(48,n,m)*(zc(kc1)-sz)**2.0d0)
        pdzz(2,n,m)=dz2*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc)-sz)**2.0d0)
      enddo

      do n=0,nickc(m)
        signy=sign(1.0,fickc(11,n,m))
        sy=fickc(3,n,m)
        je=ickc(3,n,m)
        jc=ickc(4,n,m)
        jc1=jc+1
        if(je.eq.jc) then
          pdy(n,m)=-dy1*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc1)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc1)-sy)**2.0d0)
        else
          pdy(n,m)=-dy1*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc)-sy)**2.0d0)
        endif
        pdyy(1,n,m)=-dy2*(fickc(49,n,m)*signy+
     .                    fickc(23,n,m)*(yc(jc1)-sy)+
     .              0.5d0*fickc(46,n,m)*(yc(jc1)-sy)**2.0d0)
        pdyy(2,n,m)=dy2*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc)-sy)**2.0d0)
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        sx=fjckc(3,n,m)
        ie=jckc(3,n,m)
        ic=jckc(4,n,m)
        ic1=ic+1
        if(ie.eq.ic) then
          pdx(n,m)=-dx1*(fjckc(49,n,m)*signx+
     .                   fjckc(22,n,m)*(xc(ic1)-sx)+
     .             0.5d0*fjckc(43,n,m)*(xc(ic1)-sx)**2.0d0)
        else
          pdx(n,m)=-dx1*(fjckc(49,n,m)*signx+
     .                   fjckc(22,n,m)*(xc(ic)-sx)+
     .             0.5d0*fjckc(43,n,m)*(xc(ic)-sx)**2.0d0)
        endif
        pdxx(1,n,m)=-dx2*(fjckc(49,n,m)*signx+
     .                    fjckc(22,n,m)*(xc(ic1)-sx)+
     .              0.5d0*fjckc(43,n,m)*(xc(ic1)-sx)**2.0d0)
        pdxx(2,n,m)=dx2*(fjckc(49,n,m)*signx+
     .                   fjckc(22,n,m)*(xc(ic)-sx)+
     .             0.5d0*fjckc(43,n,m)*(xc(ic)-sx)**2.0d0)
      enddo

      ENDDO

      return
      end      


c-----------------------------------------------------------------------


