c-----------------------------------------------------------------------
c
      subroutine udu_surface
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc,ie1,je1,ke1,ic1,jc1,kc1
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
          fiejc(56,n,m)=(hp*fece(ie,jc,ke)+hm*fecc(ie,jc,kc1))/hdz-
     .                fiejc(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                fiejc(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fece(ie,jc,ke1)-fece(ie,jc,ke))-
     .       dz1*(fiejc(15,n,m)*(ze(ke)-sz)+0.5d0*
     .            fiejc(30,n,m)*(ze(ke)-sz)**2.0d0)
          dm=dz1*(fecc(ie,jc,kc1)-fecc(ie,jc,kc))-
     .       dz1*(fiejc(15,n,m)*(zc(kc1)-sz)+0.5d0*
     .            fiejc(30,n,m)*(zc(kc1)-sz)**2.0d0)
          fiejc(59,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*fiejc(15,n,m)/hdz-hp*hm*fiejc(30,n,m)/hdz
          fiejc(60,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiejc(15,n,m)/hdz-hp*hm*fiejc(30,n,m)/hdz
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          fiejc(56,n,m)=(hp*fecc(ie,jc,kc)+hm*fece(ie,jc,ke1))/hdz-
     .                fiejc(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                fiejc(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
          dp=dz1*(fecc(ie,jc,kc1)-fecc(ie,jc,kc))-
     .       dz1*(fiejc(15,n,m)*(zc(kc)-sz)+0.5d0*
     .            fiejc(30,n,m)*(zc(kc)-sz)**2.0d0)
          dm=dz1*(fece(ie,jc,ke1)-fece(ie,jc,ke))-
     .       dz1*(fiejc(15,n,m)*(ze(ke1)-sz)+0.5d0*
     .            fiejc(30,n,m)*(ze(ke1)-sz)**2.0d0)
          fiejc(59,n,m)=(hp*dm+hm*dp)/hdz+
     .                  hp*fiejc(15,n,m)/hdz-hp*hm*fiejc(30,n,m)/hdz
          fiejc(60,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiejc(15,n,m)/hdz-hp*hm*fiejc(30,n,m)/hdz
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
          ficje(56,n,m)=(hp*fcee(ic,je,ke)+hm*fcec(ic,je,kc1))/hdz-
     .                ficje(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficje(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          ficje(56,n,m)=(hp*fcec(ic,je,kc)+hm*fcee(ic,je,ke1))/hdz-
     .                ficje(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficje(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
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
          ficjc(56,n,m)=(hp*fcce(ic,jc,ke)+hm*fccc(ic,jc,kc1))/hdz-
     .                ficjc(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficjc(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
        else
          hp=ze(ke1)-sz
          hm=hdz-hp
          ficjc(56,n,m)=(hp*fccc(ic,jc,kc)+hm*fcce(ic,jc,ke1))/hdz-
     .                ficjc(15,n,m)*hp*hm/hdz-signz*dz1*hp*hm*
     .                ficjc(30,n,m)*((1.0d0+signz)*hp+(1.0d0-signz)*hm)
        endif
        hp=zc(kc1)-sz
        hm=dz-hp
        dp=data1(ic,jc,kc1)
        dm=data1(ic,jc,kc)
        ficjc(59,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*signz*signx*ficjc(13,n,m)*dz1-
     .                 hp*hm*signz*ficjc(27,n,m)*dz1
        ficjc(60,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*signz*signx*ficjc(13,n,m)*dz1-
     .                 hp*hm*signz*ficjc(27,n,m)*dz1
        dp=data2(ic,jc,kc1)
        dm=data2(ic,jc,kc)
        ficjc(61,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*signz*signy*ficjc(14,n,m)*dz1-
     .                 hp*hm*signz*ficjc(29,n,m)*dz1
        ficjc(62,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*signz*signy*ficjc(14,n,m)*dz1-
     .                 hp*hm*signz*ficjc(29,n,m)*dz1
        dp=data3(ic,jc,kc1)
        dm=data3(ic,jc,kc)
        ficjc(63,n,m)=(hp*dm+hm*dp)*dz1+
     .                 hp*ficjc(15,n,m)*dz1-hp*hm*ficjc(30,n,m)*dz1
        ficjc(64,n,m)=(hp*dm+hm*dp)*dz1-
     .                 hm*ficjc(15,n,m)*dz1-hp*hm*ficjc(30,n,m)*dz1
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
          fiekc(56,n,m)=(hp*feec(ie,je,kc)+hm*fecc(ie,jc1,kc))/hdy-
     .                fiekc(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fiekc(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          dp=dy1*(feec(ie,je1,kc)-feec(ie,je,kc))-
     .       dy1*(fiekc(14,n,m)*(ye(je)-sy)+0.5d0*
     .            fiekc(28,n,m)*(ye(je)-sy)**2.0d0)
          dm=dy1*(fecc(ie,jc1,kc)-fecc(ie,jc,kc))-
     .       dy1*(fiekc(14,n,m)*(yc(jc1)-sy)+0.5d0*
     .            fiekc(28,n,m)*(yc(jc1)-sy)**2.0d0)
          fiekc(59,n,m)=(hp*dm+hm*dp)/hdy+
     .                  hp*fiekc(14,n,m)/hdy-hp*hm*fiekc(28,n,m)/hdy
          fiekc(60,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiekc(14,n,m)/hdy-hp*hm*fiekc(28,n,m)/hdy
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          fiekc(56,n,m)=(hp*fecc(ie,jc,kc)+hm*feec(ie,je1,kc))/hdy-
     .                fiekc(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fiekc(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
          dp=dy1*(fecc(ie,jc1,kc)-fecc(ie,jc,kc))-
     .       dy1*(fiekc(14,n,m)*(yc(jc)-sy)+0.5d0*
     .            fiekc(28,n,m)*(yc(jc)-sy)**2.0d0)
          dm=dy1*(feec(ie,je1,kc)-feec(ie,je,kc))-
     .       dy1*(fiekc(14,n,m)*(ye(je1)-sy)+0.5d0*
     .            fiekc(28,n,m)*(ye(je1)-sy)**2.0d0)
          fiekc(59,n,m)=(hp*dm+hm*dp)/hdy+
     .                  hp*fiekc(14,n,m)/hdy-hp*hm*fiekc(28,n,m)/hdy
          fiekc(60,n,m)=(hp*dm+hm*dp)/hdz-
     .                  hm*fiekc(14,n,m)/hdy-hp*hm*fiekc(28,n,m)/hdy
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
          ficke(56,n,m)=(hp*fcee(ic,je,ke)+hm*fcce(ic,jc1,ke))/hdy-
     .                ficke(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                ficke(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          ficke(56,n,m)=(hp*fcce(ic,jc,ke)+hm*fcee(ic,je1,ke))/hdy-
     .                ficke(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                ficke(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
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
          fickc(56,n,m)=(hp*fcec(ic,je,kc)+hm*fccc(ic,jc1,kc))/hdy-
     .                fickc(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fickc(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        else
          hp=ye(je1)-sy
          hm=hdy-hp
          fickc(56,n,m)=(hp*fccc(ic,jc,kc)+hm*fcec(ic,je1,kc))/hdy-
     .                fickc(14,n,m)*hp*hm/hdy-signy*dy1*hp*hm*
     .                fickc(28,n,m)*((1.0d0+signy)*hp+(1.0d0-signy)*hm)
        endif
        hp=yc(jc1)-sy
        hm=dy-hp
        dp=data1(ic,jc1,kc)
        dm=data1(ic,jc,kc)
        fickc(59,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*signy*signx*fickc(13,n,m)*dy1-
     .                 hp*hm*signy*fickc(26,n,m)*dy1
        fickc(60,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*signy*signx*fickc(13,n,m)*dy1-
     .                 hp*hm*signy*fickc(26,n,m)*dy1
        dp=data2(ic,jc1,kc)
        dm=data2(ic,jc,kc)
        fickc(61,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*fickc(14,n,m)*dy1-hp*hm*fickc(28,n,m)*dy1
        fickc(62,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*fickc(14,n,m)*dy1-hp*hm*fickc(28,n,m)*dy1
        dp=data3(ic,jc1,kc)
        dm=data3(ic,jc,kc)
        fickc(63,n,m)=(hp*dm+hm*dp)*dy1+
     .                 hp*signy*signz*fickc(15,n,m)*dy1-
     .                 hp*hm*signy*fickc(29,n,m)*dy1
        fickc(64,n,m)=(hp*dm+hm*dp)*dy1-
     .                 hm*signy*signz*fickc(15,n,m)*dy1-
     .                 hp*hm*signy*fickc(29,n,m)*dy1
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
          fjekc(56,n,m)=(hp*feec(ie,je,kc)+hm*fcec(ic1,je,kc))/hdx-
     .                fjekc(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjekc(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(feec(ie1,je,kc)-feec(ie,je,kc))-
     .       dx1*(fjekc(13,n,m)*(xe(ie)-sx)+0.5d0*
     .            fjekc(25,n,m)*(xe(ie)-sx)**2.0d0)
          dm=dx1*(fcec(ic1,je,kc)-fcec(ic,je,kc))-
     .       dx1*(fjekc(13,n,m)*(xc(ic1)-sx)+0.5d0*
     .            fjekc(25,n,m)*(xc(ic1)-sx)**2.0d0)
          fjekc(59,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjekc(13,n,m)/hdx-hp*hm*fjekc(25,n,m)/hdx
          fjekc(60,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjekc(13,n,m)/hdx-hp*hm*fjekc(25,n,m)/hdx
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjekc(56,n,m)=(hp*fcec(ic,je,kc)+hm*feec(ie1,je,kc))/hdx-
     .                fjekc(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjekc(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(fcec(ic1,je,kc)-fcec(ic,je,kc))-
     .       dx1*(fjekc(13,n,m)*(xc(ic)-sx)+0.5d0*
     .            fjekc(25,n,m)*(xc(ic)-sx)**2.0d0)
          dm=dx1*(feec(ie1,je,kc)-feec(ie,je,kc))-
     .       dx1*(fjekc(13,n,m)*(xe(ie1)-sx)+0.5d0*
     .            fjekc(25,n,m)*(xe(ie1)-sx)**2.0d0)
          fjekc(59,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjekc(13,n,m)/hdx-hp*hm*fjekc(25,n,m)/hdx
          fjekc(60,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjekc(13,n,m)/hdx-hp*hm*fjekc(25,n,m)/hdx
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
          fjcke(56,n,m)=(hp*fece(ie,jc,ke)+hm*fcce(ic1,jc,ke))/hdx-
     .                fjcke(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjcke(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(fece(ie1,jc,ke)-fece(ie,jc,ke))-
     .       dx1*(fjcke(13,n,m)*(xe(ie)-sx)+0.5d0*
     .            fjcke(25,n,m)*(xe(ie)-sx)**2.0d0)
          dm=dx1*(fcce(ic1,jc,ke)-fcce(ic,jc,ke))-
     .       dx1*(fjcke(13,n,m)*(xc(ic1)-sx)+0.5d0*
     .            fjcke(25,n,m)*(xc(ic1)-sx)**2.0d0)
          fjcke(59,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjcke(13,n,m)/hdx-hp*hm*fjcke(25,n,m)/hdx
          fjcke(60,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjcke(13,n,m)/hdx-hp*hm*fjcke(25,n,m)/hdx
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjcke(56,n,m)=(hp*fcce(ic,jc,ke)+hm*fece(ie1,jc,ke))/hdx-
     .                fjcke(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjcke(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
          dp=dx1*(fcce(ic1,jc,ke)-fcce(ic,jc,ke))-
     .       dx1*(fjcke(13,n,m)*(xc(ic)-sx)+0.5d0*
     .            fjcke(25,n,m)*(xc(ic)-sx)**2.0d0)
          dm=dx1*(fece(ie1,jc,ke)-fece(ie,jc,ke))-
     .       dx1*(fjcke(13,n,m)*(xe(ie1)-sx)+0.5d0*
     .            fjcke(25,n,m)*(xe(ie1)-sx)**2.0d0)
          fjcke(59,n,m)=(hp*dm+hm*dp)/hdx+
     .                  hp*fjcke(13,n,m)/hdx-hp*hm*fjcke(25,n,m)/hdx
          fjcke(60,n,m)=(hp*dm+hm*dp)/hdx-
     .                  hm*fjcke(13,n,m)/hdx-hp*hm*fjcke(25,n,m)/hdx
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
          fjckc(56,n,m)=(hp*fecc(ie,jc,kc)+hm*fccc(ic1,jc,kc))/hdx-
     .                fjckc(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjckc(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        else
          hp=xe(ie1)-sx
          hm=hdx-hp
          fjckc(56,n,m)=(hp*fccc(ic,jc,kc)+hm*fecc(ie1,jc,kc))/hdx-
     .                fjckc(13,n,m)*hp*hm/hdx-signx*dx1*hp*hm*
     .                fjckc(25,n,m)*((1.0d0+signx)*hp+(1.0d0-signx)*hm)
        endif
        hp=xc(ic1)-sx
        hm=dx-hp
        dp=data1(ic1,jc,kc)
        dm=data1(ic,jc,kc)
        fjckc(59,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*fjckc(13,n,m)*dx1-hp*hm*fjckc(25,n,m)*dx1
        fjckc(60,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*fjckc(13,n,m)*dx1-hp*hm*fjckc(25,n,m)*dx1
        dp=data2(ic1,jc,kc)
        dm=data2(ic,jc,kc)
        fjckc(61,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*signx*signy*fjckc(14,n,m)*dx1-
     .                 hp*hm*signx*fjckc(26,n,m)*dx1
        fjckc(62,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*signx*signy*fjckc(14,n,m)*dx1-
     .                 hp*hm*signx*fjckc(26,n,m)*dx1
        dp=data3(ic1,jc,kc)
        dm=data3(ic,jc,kc)
        fjckc(63,n,m)=(hp*dm+hm*dp)*dx1+
     .                 hp*signx*signz*fjckc(15,n,m)*dx1-
     .                 hp*hm*signx*fjckc(27,n,m)*dx1
        fjckc(64,n,m)=(hp*dm+hm*dp)*dx1-
     .                 hm*signx*signz*fjckc(15,n,m)*dx1-
     .                 hp*hm*signx*fjckc(27,n,m)*dx1
      enddo

      ENDDO

      return
      end


c-----------------------------------------------------------------------



