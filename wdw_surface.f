c-----------------------------------------------------------------------
c
      subroutine wdw_surface
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc,ie1,je1,ke1,ic1,jc1,kc1,np
      real*8 sx,sy,sz,hp,hm,signx,signy,signz,dp,dm

      DO m=1,ms

      do n=0,niejc(m)
        signz=sign(1.0,fiejc(12,n,m))
        sz=fiejc(3,n,m)
        ie=iejc(1,n,m)
        jc=iejc(2,n,m)
        ke=iejc(3,n,m)
        kc=iejc(4,n,m)
        ke1=ke+1
        kc1=kc+1
        if(ke.eq.kc) then
          hp=zc(kc1)-sz
          hm=hdz-hp
          fiejc(58,n,m)=(hp*fece(ie,jc,ke)+hm*fecc(ie,jc,kc1))/hdz-
     .                fiejc(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                fiejc(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fece(ie,jc,ke1)-fece(ie,jc,ke))-
     .       dz1*(fiejc(21,n,m)*(ze(ke)-sz)+0.5d0*
     .            fiejc(42,n,m)*(ze(ke)-sz)**2.0d0)
          dm=dz1*(fecc(ie,jc,kc1)-fecc(ie,jc,kc))-
     .       dz1*(fiejc(21,n,m)*(zc(kc1)-sz)+0.5d0*
     .            fiejc(42,n,m)*(zc(kc1)-sz)**2.0d0)
          fiejc(61,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*fiejc(21,n,m)/hdz-hp*hm*fiejc(42,n,m)/hdz
          fiejc(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiejc(21,n,m)/hdz-hp*hm*fiejc(42,n,m)/hdz
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          fiejc(58,n,m)=(hp*fecc(ie,jc,kc)+hm*fece(ie,jc,ke1))/hdz-
     .                fiejc(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                fiejc(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fecc(ie,jc,kc1)-fecc(ie,jc,kc))-
     .       dz1*(fiejc(21,n,m)*(zc(kc)-sz)+0.5d0*
     .            fiejc(42,n,m)*(zc(kc)-sz)**2.0d0)
          dm=dz1*(fece(ie,jc,ke1)-fece(ie,jc,ke))-
     .       dz1*(fiejc(21,n,m)*(ze(ke1)-sz)+0.5d0*
     .            fiejc(42,n,m)*(ze(ke1)-sz)**2.0d0)
          fiejc(61,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*fiejc(21,n,m)/hdz-hp*hm*fiejc(42,n,m)/hdz
          fiejc(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiejc(21,n,m)/hdz-hp*hm*fiejc(42,n,m)/hdz
        endif
      enddo

      do n=0,nicje(m)
        signz=sign(1.0,ficje(12,n,m))
        sz=ficje(3,n,m)
        ic=icje(1,n,m)
        je=icje(2,n,m)
        ke=icje(3,n,m)
        kc=icje(4,n,m)
        ke1=ke+1
        kc1=kc+1
        if(ke.eq.kc) then
          hp=zc(kc1)-sz
          hm=hdz-hp
          ficje(58,n,m)=(hp*fcee(ic,je,ke)+hm*fcec(ic,je,kc1))/hdz-
     .                ficje(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficje(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fcee(ic,je,ke1)-fcee(ic,je,ke))-
     .       dz1*(ficje(21,n,m)*(ze(ke)-sz)+0.5d0*
     .            ficje(42,n,m)*(ze(ke)-sz)**2.0d0)
          dm=dz1*(fcec(ic,je,kc1)-fcec(ic,je,kc))-
     .       dz1*(ficje(21,n,m)*(zc(kc1)-sz)+0.5d0*
     .            ficje(42,n,m)*(zc(kc1)-sz)**2.0d0)
          ficje(61,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*ficje(21,n,m)/hdz-hp*hm*ficje(42,n,m)/hdz
          ficje(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*ficje(21,n,m)/hdz-hp*hm*ficje(42,n,m)/hdz
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          ficje(58,n,m)=(hp*fcec(ic,je,kc)+hm*fcee(ic,je,ke1))/hdz-
     .                ficje(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficje(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fcec(ic,je,kc1)-fcec(ic,je,kc))-
     .       dz1*(ficje(21,n,m)*(zc(kc)-sz)+0.5d0*
     .            ficje(42,n,m)*(zc(kc)-sz)**2.0d0)
          dm=dz1*(fcee(ic,je,ke1)-fcee(ic,je,ke))-
     .       dz1*(ficje(21,n,m)*(ze(ke1)-sz)+0.5d0*
     .            ficje(42,n,m)*(ze(ke1)-sz)**2.0d0)
          ficje(61,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*ficje(21,n,m)/hdz-hp*hm*ficje(42,n,m)/hdz
          ficje(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*ficje(21,n,m)/hdz-hp*hm*ficje(42,n,m)/hdz
        endif
      enddo

      do n=0,nicjc(m)
        signx=sign(1.0,ficjc(10,n,m))
        signy=sign(1.0,ficjc(11,n,m))
        signz=sign(1.0,ficjc(12,n,m))
        sz=ficjc(3,n,m)
        ic=icjc(1,n,m)
        jc=icjc(2,n,m)
        ke=icjc(3,n,m)
        kc=icjc(4,n,m)
        ke1=ke+1
        kc1=kc+1
        if(ke.eq.kc) then
          hp=zc(kc1)-sz
          hm=hdz-hp
          ficjc(58,n,m)=(hp*fcce(ic,jc,ke)+hm*fccc(ic,jc,kc1))/hdz-
     .                ficjc(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficjc(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          ficjc(58,n,m)=(hp*fccc(ic,jc,kc)+hm*fcce(ic,jc,ke1))/hdz-
     .                ficjc(21,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficjc(42,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
        endif
        hp=zc(kc1)-sz
        hm=dz-hp
        dp=data7(ic,jc,kc1)
        dm=data7(ic,jc,kc)
        ficjc(71,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*signz*signx*ficjc(19,n,m)*dz1-
     .                 hp*hm*signz*ficjc(39,n,m)*dz1
        ficjc(72,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*signz*signx*ficjc(19,n,m)*dz1-
     .                 hp*hm*signz*ficjc(39,n,m)*dz1
        dp=data8(ic,jc,kc1)
        dm=data8(ic,jc,kc)
        ficjc(73,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*signz*signy*ficjc(20,n,m)*dz1-
     .                 hp*hm*signz*ficjc(41,n,m)*dz1
        ficjc(74,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*signz*signy*ficjc(20,n,m)*dz1-
     .                 hp*hm*signz*ficjc(41,n,m)*dz1
        dp=data9(ic,jc,kc1)
        dm=data9(ic,jc,kc)
        ficjc(75,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*ficjc(21,n,m)*dz1-hp*hm*ficjc(42,n,m)*dz1
        ficjc(76,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*ficjc(21,n,m)*dz1-hp*hm*ficjc(42,n,m)*dz1
      enddo

      do n=0,niekc(m)
        signy=sign(1.0,fiekc(11,n,m))
        sy=fiekc(3,n,m)
        ie=iekc(1,n,m)
        kc=iekc(2,n,m)
        je=iekc(3,n,m)
        jc=iekc(4,n,m)
        je1=je+1
        jc1=jc+1
        if(je.eq.jc) then
          hp=yc(jc1)-sy
          hm=hdy-hp
          fiekc(58,n,m)=(hp*feec(ie,je,kc)+hm*fecc(ie,jc1,kc))/hdy-
     .                fiekc(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fiekc(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          fiekc(58,n,m)=(hp*fecc(ie,jc,kc)+hm*feec(ie,je1,kc))/hdy-
     .                fiekc(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fiekc(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        endif
      enddo

      do n=0,nicke(m)
        signy=sign(1.0,ficke(11,n,m))
        sy=ficke(3,n,m)
        ic=icke(1,n,m)
        ke=icke(2,n,m)
        je=icke(3,n,m)
        jc=icke(4,n,m)
        je1=je+1
        jc1=jc+1
        if(je.eq.jc) then
          hp=yc(jc1)-sy
          hm=hdy-hp
          ficke(58,n,m)=(hp*fcee(ic,je,ke)+hm*fcce(ic,jc1,ke))/hdy-
     .                ficke(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                ficke(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          dp=dy1*(fcee(ic,je1,ke)-fcee(ic,je,ke))-
     .       dy1*(ficke(20,n,m)*(ye(je)-sy)+0.5d0*
     .            ficke(40,n,m)*(ye(je)-sy)**2.0d0)
          dm=dy1*(fcce(ic,jc1,ke)-fcce(ic,jc,ke))-
     .       dy1*(ficke(20,n,m)*(yc(jc1)-sy)+0.5d0*
     .            ficke(40,n,m)*(yc(jc1)-sy)**2.0d0)
          ficke(61,n,m)=(hp*dm+hm*dp)/hdy+
     .                  hp*ficke(20,n,m)/hdy-hp*hm*ficke(40,n,m)/hdy
          ficke(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*ficke(20,n,m)/hdy-hp*hm*ficke(40,n,m)/hdy
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          ficke(58,n,m)=(hp*fcce(ic,jc,ke)+hm*fcee(ic,je1,ke))/hdy-
     .                ficke(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                ficke(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          dp=dy1*(fcce(ic,jc1,ke)-fcce(ic,jc,ke))-
     .       dy1*(ficke(20,n,m)*(yc(jc)-sy)+0.5d0*
     .            ficke(40,n,m)*(yc(jc)-sy)**2.0d0)
          dm=dy1*(fcee(ic,je1,ke)-fcee(ic,je,ke))-
     .       dy1*(ficke(20,n,m)*(ye(je1)-sy)+0.5d0*
     .            ficke(40,n,m)*(ye(je1)-sy)**2.0d0)
          ficke(61,n,m)=(hp*dm+hm*dp)/hdy+
     .                  hp*ficke(20,n,m)/hdy-hp*hm*ficke(40,n,m)/hdy
          ficke(62,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*ficke(20,n,m)/hdy-hp*hm*ficke(40,n,m)/hdy
        endif
      enddo

      do n=0,nickc(m)
        signy=sign(1.0,fickc(11,n,m))
        sy=fickc(3,n,m)
        ic=ickc(1,n,m)
        kc=ickc(2,n,m)
        je=ickc(3,n,m)
        jc=ickc(4,n,m)
        je1=je+1
        jc1=jc+1
        if(je.eq.jc) then
          hp=yc(jc1)-sy
          hm=hdy-hp
          fickc(58,n,m)=(hp*fcec(ic,je,kc)+hm*fccc(ic,jc1,kc))/hdy-
     .                fickc(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fickc(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          fickc(58,n,m)=(hp*fccc(ic,jc,kc)+hm*fcec(ic,je1,kc))/hdy-
     .                fickc(20,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fickc(40,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        endif
        hp=yc(jc1)-sy
        hm=dy-hp
        dp=data7(ic,jc1,kc)
        dm=data7(ic,jc,kc)
        fickc(71,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*signy*signx*fickc(19,n,m)*dy1-
     .                 hp*hm*signy*fickc(38,n,m)*dy1
        fickc(72,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*signy*signx*fickc(19,n,m)*dy1-
     .                 hp*hm*signy*fickc(38,n,m)*dy1
        dp=data8(ic,jc1,kc)
        dm=data8(ic,jc,kc)
        fickc(73,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*fickc(20,n,m)*dy1-hp*hm*fickc(40,n,m)*dy1
        fickc(74,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*fickc(20,n,m)*dy1-hp*hm*fickc(40,n,m)*dy1
        dp=data9(ic,jc1,kc)
        dm=data9(ic,jc,kc)
        fickc(75,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*signy*signz*fickc(21,n,m)*dy1-
     .                 hp*hm*signy*fickc(41,n,m)*dy1
        fickc(76,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*signy*signz*fickc(21,n,m)*dy1-
     .                 hp*hm*signy*fickc(41,n,m)*dy1
      enddo

      do n=0,njekc(m)
        signx=sign(1.0,fjekc(10,n,m))
        sx=fjekc(3,n,m)
        je=jekc(1,n,m)
        kc=jekc(2,n,m)
        ie=jekc(3,n,m)
        ic=jekc(4,n,m)
        ie1=ie+1
        ic1=ic+1
        if(ie.eq.ic) then
          hp=xc(ic1)-sx
          hm=hdx-hp
          fjekc(58,n,m)=(hp*feec(ie,je,kc)+hm*fcec(ic1,je,kc))/hdx-
     .                fjekc(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjekc(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjekc(58,n,m)=(hp*fcec(ic,je,kc)+hm*feec(ie1,je,kc))/hdx-
     .                fjekc(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjekc(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        endif
      enddo

      do n=0,njcke(m)
        signx=sign(1.0,fjcke(10,n,m))
        sx=fjcke(3,n,m)
        jc=jcke(1,n,m)
        ke=jcke(2,n,m)
        ie=jcke(3,n,m)
        ic=jcke(4,n,m)
        ie1=ie+1
        ic1=ic+1
        if(ie.eq.ic) then
          hp=xc(ic1)-sx
          hm=hdx-hp
          fjcke(58,n,m)=(hp*fece(ie,jc,ke)+hm*fcce(ic1,jc,ke))/hdx-
     .                fjcke(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjcke(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(fece(ie1,jc,ke)-fece(ie,jc,ke))-
     .       dx1*(fjcke(19,n,m)*(xe(ie)-sx)+0.5d0*
     .            fjcke(37,n,m)*(xe(ie)-sx)**2.0d0)
          dm=dx1*(fcce(ic1,jc,ke)-fcce(ic,jc,ke))-
     .       dx1*(fjcke(19,n,m)*(xc(ic1)-sx)+0.5d0*
     .            fjcke(37,n,m)*(xc(ic1)-sx)**2.0d0)
          fjcke(61,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjcke(19,n,m)/hdx-hp*hm*fjcke(37,n,m)/hdx
          fjcke(62,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjcke(19,n,m)/hdx-hp*hm*fjcke(37,n,m)/hdx
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjcke(58,n,m)=(hp*fcce(ic,jc,ke)+hm*fece(ie1,jc,ke))/hdx-
     .                fjcke(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjcke(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(fcce(ic1,jc,ke)-fcce(ic,jc,ke))-
     .       dx1*(fjcke(19,n,m)*(xc(ic)-sx)+0.5d0*
     .            fjcke(37,n,m)*(xc(ic)-sx)**2.0d0)
          dm=dx1*(fece(ie1,jc,ke)-fece(ie,jc,ke))-
     .       dx1*(fjcke(19,n,m)*(xe(ie1)-sx)+0.5d0*
     .            fjcke(37,n,m)*(xe(ie1)-sx)**2.0d0)
          fjcke(61,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjcke(19,n,m)/hdx-hp*hm*fjcke(37,n,m)/hdx
          fjcke(62,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjcke(19,n,m)/hdx-hp*hm*fjcke(37,n,m)/hdx
        endif
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        sx=fjckc(3,n,m)
        jc=jckc(1,n,m)
        kc=jckc(2,n,m)
        ie=jckc(3,n,m)
        ic=jckc(4,n,m)
        ie1=ie+1
        ic1=ic+1
        if(ie.eq.ic) then
          hp=xc(ic1)-sx
          hm=hdx-hp
          fjckc(58,n,m)=(hp*fecc(ie,jc,kc)+hm*fccc(ic1,jc,kc))/hdx-
     .                fjckc(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjckc(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjckc(58,n,m)=(hp*fccc(ic,jc,kc)+hm*fecc(ie1,jc,kc))/hdx-
     .                fjckc(19,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjckc(37,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        endif
        hp=xc(ic1)-sx
        hm=dx-hp
        dp=data7(ic1,jc,kc)
        dm=data7(ic,jc,kc)
        fjckc(71,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*fjckc(19,n,m)*dx1-hp*hm*fjckc(37,n,m)*dx1
        fjckc(72,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*fjckc(19,n,m)*dx1-hp*hm*fjckc(37,n,m)*dx1
        dp=data8(ic1,jc,kc)
        dm=data8(ic,jc,kc)
        fjckc(73,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*signx*signy*fjckc(20,n,m)*dx1-
     .                 hp*hm*signx*fjckc(38,n,m)*dx1
        fjckc(74,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*signx*signy*fjckc(20,n,m)*dx1-
     .                 hp*hm*signx*fjckc(38,n,m)*dx1
        dp=data9(ic1,jc,kc)
        dm=data9(ic,jc,kc)
        fjckc(75,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*signx*signz*fjckc(21,n,m)*dx1-
     .                 hp*hm*signx*fjckc(39,n,m)*dx1
        fjckc(76,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*signx*signz*fjckc(21,n,m)*dx1-
     .                 hp*hm*signx*fjckc(39,n,m)*dx1
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------



