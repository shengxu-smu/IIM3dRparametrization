c-----------------------------------------------------------------------
c
      subroutine correction_difference
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc,ie1,je1,ke1,ic1,jc1,kc1
      real*8 sx,sy,sz,su,sv,sw,djc,ddjc,signx,signy,signz

      DO m=1,ms

      do n=0,niejc(m)
        sz=fiejc(3,n,m)
        su=fiejc(56,n,m)
        sv=fiejc(57,n,m)
        sw=fiejc(58,n,m)
        ke=iejc(3,n,m)
        ke1=ke+1
        kc=iejc(4,n,m)
        kc1=kc+1
        djc=su*fiejc(21,n,m)+sw*fiejc(15,n,m)
        ddjc=su*fiejc(42,n,m)+sw*fiejc(30,n,m)+2.0d0*
     .       (fiejc(59,n,m)*fiejc(61,n,m)-fiejc(60,n,m)*fiejc(62,n,m))
        if(ke.eq.kc) then
          uwdz(n,m)=-dz1*(djc*(ze(ke)-sz)+
     .             0.5d0*ddjc*(ze(ke)-sz)**2.0d0)
        else
          uwdz(n,m)=-dz1*(djc*(ze(ke1)-sz)+
     .             0.5d0*ddjc*(ze(ke1)-sz)**2.0d0)
        endif
        udzz(1,n,m)=-dz2*(fiejc(15,n,m)*(zc(kc1)-sz)+
     .              0.5d0*fiejc(30,n,m)*(zc(kc1)-sz)**2.0d0)
        udzz(2,n,m)=dz2*(fiejc(15,n,m)*(zc(kc)-sz)+
     .             0.5d0*fiejc(30,n,m)*(zc(kc)-sz)**2.0d0)
      enddo

      do n=0,nicje(m)
        sz=ficje(3,n,m)
        su=ficje(56,n,m)
        sv=ficje(57,n,m)
        sw=ficje(58,n,m)
        ke=icje(3,n,m)
        ke1=ke+1
        kc=icje(4,n,m)
        kc1=kc+1
        djc=sv*ficje(21,n,m)+sw*ficje(18,n,m)
        ddjc=sv*ficje(42,n,m)+sw*ficje(36,n,m)+2.0d0*
     .       (ficje(59,n,m)*ficje(61,n,m)-ficje(60,n,m)*ficje(62,n,m))
        if(ke.eq.kc) then
          vwdz(n,m)=-dz1*(djc*(ze(ke)-sz)+
     .             0.5d0*ddjc*(ze(ke)-sz)**2.0d0)
        else
          vwdz(n,m)=-dz1*(djc*(ze(ke1)-sz)+
     .             0.5d0*ddjc*(ze(ke1)-sz)**2.0d0)
        endif
        vdzz(1,n,m)=-dz2*(ficje(18,n,m)*(zc(kc1)-sz)+
     .              0.5d0*ficje(36,n,m)*(zc(kc1)-sz)**2.0d0)
        vdzz(2,n,m)=dz2*(ficje(18,n,m)*(zc(kc)-sz)+
     .             0.5d0*ficje(36,n,m)*(zc(kc)-sz)**2.0d0)
      enddo

      do n=0,nicjc(m)
        signz=sign(1.0,ficjc(12,n,m))
        sz=ficjc(3,n,m)
        su=ficjc(56,n,m)
        sv=ficjc(57,n,m)
        sw=ficjc(58,n,m)
        ke=icjc(3,n,m)
        ke1=ke+1
        kc=icjc(4,n,m)
        kc1=kc+1
        djc=2.0d0*sw*ficjc(21,n,m)
        ddjc=2.0d0*sw*ficjc(42,n,m)+2.0d0*
     .       (ficjc(75,n,m)*ficjc(75,n,m)-ficjc(76,n,m)*ficjc(76,n,m))
        if(ke.eq.kc) then
          wwdz(n,m)=-dz1*(djc*(zc(kc1)-sz)+
     .             0.5d0*ddjc*(zc(kc1)-sz)**2.0d0)
          pdz(n,m)=-dz1*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc1)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc1)-sz)**2.0d0)
        else
          wwdz(n,m)=-dz1*(djc*(zc(kc)-sz)+
     .             0.5d0*ddjc*(zc(kc)-sz)**2.0d0)
          pdz(n,m)=-dz1*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc)-sz)**2.0d0)
        endif
        wdzz(1,n,m)=-dz2*(ficjc(21,n,m)*(ze(ke1)-sz)+
     .              0.5d0*ficjc(42,n,m)*(ze(ke1)-sz)**2.0d0)
        wdzz(2,n,m)=dz2*(ficjc(21,n,m)*(ze(ke)-sz)+
     .             0.5d0*ficjc(42,n,m)*(ze(ke)-sz)**2.0d0)
        pdzz(1,n,m)=-dz2*(ficjc(49,n,m)*signz+
     .                    ficjc(24,n,m)*(zc(kc1)-sz)+
     .              0.5d0*ficjc(48,n,m)*(zc(kc1)-sz)**2.0d0)
        pdzz(2,n,m)=dz2*(ficjc(49,n,m)*signz+
     .                   ficjc(24,n,m)*(zc(kc)-sz)+
     .             0.5d0*ficjc(48,n,m)*(zc(kc)-sz)**2.0d0)
      enddo

      do n=0,niekc(m)
        sy=fiekc(3,n,m)
        su=fiekc(56,n,m)
        sv=fiekc(57,n,m)
        sw=fiekc(58,n,m)
        je=iekc(3,n,m)
        je1=je+1
        jc=iekc(4,n,m)
        jc1=jc+1
        djc=su*fiekc(17,n,m)+sv*fiekc(14,n,m)
        ddjc=su*fiekc(34,n,m)+sv*fiekc(28,n,m)+2.0d0*
     .       (fiekc(59,n,m)*fiekc(61,n,m)-fiekc(60,n,m)*fiekc(62,n,m))
        if(je.eq.jc) then
          uvdy(n,m)=-dy1*(djc*(ye(je)-sy)+
     .             0.5d0*ddjc*(ye(je)-sy)**2.0d0)
        else
          uvdy(n,m)=-dy1*(djc*(ye(je1)-sy)+
     .             0.5d0*ddjc*(ye(je1)-sy)**2.0d0)
        endif
        udyy(1,n,m)=-dy2*(fiekc(14,n,m)*(yc(jc1)-sy)+
     .              0.5d0*fiekc(28,n,m)*(yc(jc1)-sy)**2.0d0)
        udyy(2,n,m)=dy2*(fiekc(14,n,m)*(yc(jc)-sy)+
     .             0.5d0*fiekc(28,n,m)*(yc(jc)-sy)**2.0d0)
      enddo

      do n=0,nicke(m)
        sy=ficke(3,n,m)
        su=ficke(56,n,m)
        sv=ficke(57,n,m)
        sw=ficke(58,n,m)
        je=icke(3,n,m)
        je1=je+1
        jc=icke(4,n,m)
        jc1=jc+1
        djc=sv*ficke(20,n,m)+sw*ficke(17,n,m)
        ddjc=sv*ficke(40,n,m)+sw*ficke(34,n,m)+2.0d0*
     .       (ficke(59,n,m)*ficke(61,n,m)-ficke(60,n,m)*ficke(62,n,m))
        if(je.eq.jc) then
          wvdy(n,m)=-dy1*(djc*(ye(je)-sy)+
     .             0.5d0*ddjc*(ye(je)-sy)**2.0d0)
        else
          wvdy(n,m)=-dy1*(djc*(ye(je1)-sy)+
     .             0.5d0*ddjc*(ye(je1)-sy)**2.0d0)
        endif
        wdyy(1,n,m)=-dy2*(ficke(20,n,m)*(yc(jc1)-sy)+
     .              0.5d0*ficke(40,n,m)*(yc(jc1)-sy)**2.0d0)
        wdyy(2,n,m)=dy2*(ficke(20,n,m)*(yc(jc)-sy)+
     .             0.5d0*ficke(40,n,m)*(yc(jc)-sy)**2.0d0)
      enddo

      do n=0,nickc(m)
        signy=sign(1.0,fickc(11,n,m))
        sy=fickc(3,n,m)
        su=fickc(56,n,m)
        sv=fickc(57,n,m)
        sw=fickc(58,n,m)
        je=ickc(3,n,m)
        je1=je+1
        jc=ickc(4,n,m)
        jc1=jc+1
        djc=2.0d0*sv*fickc(17,n,m)
        ddjc=2.0d0*sv*fickc(34,n,m)+2.0d0*
     .       (fickc(67,n,m)*fickc(67,n,m)-fickc(68,n,m)*fickc(68,n,m))
        if(je.eq.jc) then
          vvdy(n,m)=-dy1*(djc*(yc(jc1)-sy)+
     .             0.5d0*ddjc*(yc(jc1)-sy)**2.0d0)
          pdy(n,m)=-dy1*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc1)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc1)-sy)**2.0d0)
        else
          vvdy(n,m)=-dy1*(djc*(yc(jc)-sy)+
     .             0.5d0*ddjc*(yc(jc)-sy)**2.0d0)
          pdy(n,m)=-dy1*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc)-sy)**2.0d0)
        endif
        vdyy(1,n,m)=-dy2*(fickc(17,n,m)*(ye(je1)-sy)+
     .              0.5d0*fickc(34,n,m)*(ye(je1)-sy)**2.0d0)
        vdyy(2,n,m)=dy2*(fickc(17,n,m)*(ye(je)-sy)+
     .             0.5d0*fickc(34,n,m)*(ye(je)-sy)**2.0d0)
        pdyy(1,n,m)=-dy2*(fickc(49,n,m)*signy+
     .                    fickc(23,n,m)*(yc(jc1)-sy)+
     .              0.5d0*fickc(46,n,m)*(yc(jc1)-sy)**2.0d0)
        pdyy(2,n,m)=dy2*(fickc(49,n,m)*signy+
     .                   fickc(23,n,m)*(yc(jc)-sy)+
     .             0.5d0*fickc(46,n,m)*(yc(jc)-sy)**2.0d0)
      enddo

      do n=0,njekc(m)
        sx=fjekc(3,n,m)
        su=fjekc(56,n,m)
        sv=fjekc(57,n,m)
        sw=fjekc(58,n,m)
        ie=jekc(3,n,m)
        ie1=ie+1
        ic=jekc(4,n,m)
        ic1=ic+1
        djc=su*fjekc(16,n,m)+sv*fjekc(13,n,m)
        ddjc=su*fjekc(31,n,m)+sv*fjekc(25,n,m)+2.0d0*
     .       (fjekc(59,n,m)*fjekc(61,n,m)-fjekc(60,n,m)*fjekc(62,n,m))
        if(ie.eq.ic) then
          vudx(n,m)=-dx1*(djc*(xe(ie)-sx)+
     .             0.5d0*ddjc*(xe(ie)-sx)**2.0d0)
        else
          vudx(n,m)=-dx1*(djc*(xe(ie1)-sx)+
     .             0.5d0*ddjc*(xe(ie1)-sx)**2.0d0)
        endif
        vdxx(1,n,m)=-dx2*(fjekc(16,n,m)*(xc(ic1)-sx)+
     .              0.5d0*fjekc(31,n,m)*(xc(ic1)-sx)**2.0d0)
        vdxx(2,n,m)=dx2*(fjekc(16,n,m)*(xc(ic)-sx)+
     .             0.5d0*fjekc(31,n,m)*(xc(ic)-sx)**2.0d0)
      enddo

      do n=0,njcke(m)
        sx=fjcke(3,n,m)
        su=fjcke(56,n,m)
        sv=fjcke(57,n,m)
        sw=fjcke(58,n,m)
        ie=jcke(3,n,m)
        ie1=ie+1
        ic=jcke(4,n,m)
        ic1=ic+1
        djc=su*fjcke(19,n,m)+sw*fjcke(13,n,m)
        ddjc=su*fjcke(37,n,m)+sw*fjcke(25,n,m)+2.0d0*
     .       (fjcke(59,n,m)*fjcke(61,n,m)-fjcke(60,n,m)*fjcke(62,n,m))
        if(ie.eq.ic) then
          wudx(n,m)=-dx1*(djc*(xe(ie)-sx)+
     .             0.5d0*ddjc*(xe(ie)-sx)**2.0d0)
        else
          wudx(n,m)=-dx1*(djc*(xe(ie1)-sx)+
     .             0.5d0*ddjc*(xe(ie1)-sx)**2.0d0)
        endif
        wdxx(1,n,m)=-dx2*(fjcke(19,n,m)*(xc(ic1)-sx)+
     .              0.5d0*fjcke(37,n,m)*(xc(ic1)-sx)**2.0d0)
        wdxx(2,n,m)=dx2*(fjcke(19,n,m)*(xc(ic)-sx)+
     .             0.5d0*fjcke(37,n,m)*(xc(ic)-sx)**2.0d0)
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        sx=fjckc(3,n,m)
        su=fjckc(56,n,m)
        sv=fjckc(57,n,m)
        sw=fjckc(58,n,m)
        ie=jckc(3,n,m)
        ie1=ie+1
        ic=jckc(4,n,m)
        ic1=ic+1
        djc=2.0d0*su*fjckc(13,n,m)
        ddjc=2.0d0*su*fjckc(25,n,m)+2.0d0*
     .       (fjckc(59,n,m)*fjckc(59,n,m)-fjckc(60,n,m)*fjckc(60,n,m))
        if(ie.eq.ic) then
          uudx(n,m)=-dx1*(djc*(xc(ic1)-sx)+
     .             0.5d0*ddjc*(xc(ic1)-sx)**2.0d0)
          pdx(n,m)=-dx1*(fjckc(49,n,m)*signx+
     .                   fjckc(22,n,m)*(xc(ic1)-sx)+
     .             0.5d0*fjckc(43,n,m)*(xc(ic1)-sx)**2.0d0)
        else
          uudx(n,m)=-dx1*(djc*(xc(ic)-sx)+
     .             0.5d0*ddjc*(xc(ic)-sx)**2.0d0)
          pdx(n,m)=-dx1*(fjckc(49,n,m)*signx+
     .                   fjckc(22,n,m)*(xc(ic)-sx)+
     .             0.5d0*fjckc(43,n,m)*(xc(ic)-sx)**2.0d0)
        endif
        udxx(1,n,m)=-dx2*(fjckc(13,n,m)*(xe(ie1)-sx)+
     .              0.5d0*fjckc(25,n,m)*(xe(ie1)-sx)**2.0d0)
        udxx(2,n,m)=dx2*(fjckc(13,n,m)*(xe(ie)-sx)+
     .             0.5d0*fjckc(25,n,m)*(xe(ie)-sx)**2.0d0)
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


