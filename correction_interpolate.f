c-----------------------------------------------------------------------
c
      subroutine correction_interpolate
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer is,js,ks,ie,je,ke,ic,jc,kc
      real*8 sx,sy,sz

      DO m=1,ms

      do n=0,niejc(m)
        sz=fiejc(3,n,m)
        ke=iejc(3,n,m)
        kc=iejc(4,n,m)
        if(ke.eq.kc) then
          kuiejc(n,m)=ke
          ks=kc+1
          ukiejc(n,m)=-0.5d0*(fiejc(15,n,m)*(zc(ks)-sz)+
     .                  0.5d0*fiejc(30,n,m)*(zc(ks)-sz)**2.0d0)
          kwiejc(n,m)=kc+1
          ks=ke
          wkiejc(n,m)=0.5d0*(fiejc(21,n,m)*(ze(ks)-sz)+
     .                 0.5d0*fiejc(42,n,m)*(ze(ks)-sz)**2.0d0)          
        else
          kuiejc(n,m)=ke+1
          ks=kc
          ukiejc(n,m)=0.5d0*(fiejc(15,n,m)*(zc(ks)-sz)+
     .                 0.5d0*fiejc(30,n,m)*(zc(ks)-sz)**2.0d0)
          kwiejc(n,m)=kc
          ks=ke+1
          wkiejc(n,m)=-0.5d0*(fiejc(21,n,m)*(ze(ks)-sz)+
     .                  0.5d0*fiejc(42,n,m)*(ze(ks)-sz)**2.0d0)
        endif
      enddo

      do n=0,nicje(m)
        sz=ficje(3,n,m)
        ke=icje(3,n,m)
        kc=icje(4,n,m)
        if(ke.eq.kc) then
          kvicje(n,m)=ke
          ks=kc+1
          vkicje(n,m)=-0.5d0*(ficje(18,n,m)*(zc(ks)-sz)+
     .                  0.5d0*ficje(36,n,m)*(zc(ks)-sz)**2.0d0)
          kwicje(n,m)=kc+1
          ks=ke
          wkicje(n,m)=0.5d0*(ficje(21,n,m)*(ze(ks)-sz)+
     .                 0.5d0*ficje(42,n,m)*(ze(ks)-sz)**2.0d0)
          kuicje(n,m)=ke
          ks=kc+1
          ukicje(n,m)=-0.5d0*(ficje(15,n,m)*(zc(ks)-sz)+
     .                  0.5d0*ficje(30,n,m)*(zc(ks)-sz)**2.0d0)
        else
          kvicje(n,m)=ke+1
          ks=kc
          vkicje(n,m)=0.5d0*(ficje(18,n,m)*(zc(ks)-sz)+
     .                 0.5d0*ficje(36,n,m)*(zc(ks)-sz)**2.0d0)
          kwicje(n,m)=kc
          ks=ke+1
          wkicje(n,m)=-0.5d0*(ficje(21,n,m)*(ze(ks)-sz)+
     .                  0.5d0*ficje(42,n,m)*(ze(ks)-sz)**2.0d0)
          kuicje(n,m)=ke+1
          ks=kc
          ukicje(n,m)=0.5d0*(ficje(15,n,m)*(zc(ks)-sz)+
     .                 0.5d0*ficje(30,n,m)*(zc(ks)-sz)**2.0d0)

        endif
      enddo

      do n=0,nicjc(m)
        sz=ficjc(3,n,m)
        ke=icjc(3,n,m)
        kc=icjc(4,n,m)
        if(ke.eq.kc) then
          kwicjc(n,m)=kc+1
          ks=ke
          wkicjc(n,m)=0.5d0*(ficjc(21,n,m)*(ze(ks)-sz)+
     .                 0.5d0*ficjc(42,n,m)*(ze(ks)-sz)**2.0d0)
        else
          kwicjc(n,m)=kc
          ks=ke+1
          wkicjc(n,m)=-0.5d0*(ficjc(21,n,m)*(ze(ks)-sz)+
     .                  0.5d0*ficjc(42,n,m)*(ze(ks)-sz)**2.0d0)
        endif
      enddo

      do n=0,niekc(m)
        sy=fiekc(3,n,m)
        je=iekc(3,n,m)
        jc=iekc(4,n,m)
        if(je.eq.jc) then
          juiekc(n,m)=je
          js=jc+1
          ujiekc(n,m)=-0.5d0*(fiekc(14,n,m)*(yc(js)-sy)+
     .                  0.5d0*fiekc(28,n,m)*(yc(js)-sy)**2.0d0)
          jviekc(n,m)=jc+1
          js=je
          vjiekc(n,m)=0.5d0*(fiekc(17,n,m)*(ye(js)-sy)+
     .                 0.5d0*fiekc(34,n,m)*(ye(js)-sy)**2.0d0)
          jwiekc(n,m)=je
          js=jc+1
          wjiekc(n,m)=-0.5d0*(fiekc(20,n,m)*(yc(js)-sy)+
     .                  0.5d0*fiekc(40,n,m)*(yc(js)-sy)**2.0d0)
        else
          juiekc(n,m)=je+1
          js=jc
          ujiekc(n,m)=0.5d0*(fiekc(14,n,m)*(yc(js)-sy)+
     .                 0.5d0*fiekc(28,n,m)*(yc(js)-sy)**2.0d0)
          jviekc(n,m)=jc
          js=je+1
          vjiekc(n,m)=-0.5d0*(fiekc(17,n,m)*(ye(js)-sy)+
     .                  0.5d0*fiekc(34,n,m)*(ye(js)-sy)**2.0d0)
          jwiekc(n,m)=je+1
          js=jc
          wjiekc(n,m)=0.5d0*(fiekc(20,n,m)*(yc(js)-sy)+
     .                 0.5d0*fiekc(40,n,m)*(yc(js)-sy)**2.0d0)
        endif
      enddo

      do n=0,nicke(m)
        sy=ficke(3,n,m)
        je=icke(3,n,m)
        jc=icke(4,n,m)
        if(je.eq.jc) then
          jwicke(n,m)=je
          js=jc+1
          wjicke(n,m)=-0.5d0*(ficke(20,n,m)*(yc(js)-sy)+
     .                  0.5d0*ficke(40,n,m)*(yc(js)-sy)**2.0d0)
          jvicke(n,m)=jc+1
          js=je
          vjicke(n,m)=0.5d0*(ficke(17,n,m)*(ye(js)-sy)+
     .                 0.5d0*ficke(34,n,m)*(ye(js)-sy)**2.0d0)
        else
          jwicke(n,m)=je+1
          js=jc
          wjicke(n,m)=0.5d0*(ficke(20,n,m)*(yc(js)-sy)+
     .                 0.5d0*ficke(40,n,m)*(yc(js)-sy)**2.0d0)
          jvicke(n,m)=jc
          js=je+1
          vjicke(n,m)=-0.5d0*(ficke(17,n,m)*(ye(js)-sy)+
     .                  0.5d0*ficke(34,n,m)*(ye(js)-sy)**2.0d0)
        endif
      enddo

      do n=0,nickc(m)
        sy=fickc(3,n,m)
        je=ickc(3,n,m)
        jc=ickc(4,n,m)
        if(je.eq.jc) then
          jvickc(n,m)=jc+1
          js=je
          vjickc(n,m)=0.5d0*(fickc(17,n,m)*(ye(js)-sy)+
     .                 0.5d0*fickc(34,n,m)*(ye(js)-sy)**2.0d0)
        else
          jvickc(n,m)=jc
          js=je+1
          vjickc(n,m)=-0.5d0*(fickc(17,n,m)*(ye(js)-sy)+
     .                  0.5d0*fickc(34,n,m)*(ye(js)-sy)**2.0d0)
        endif
      enddo

      do n=0,njekc(m)
        sx=fjekc(3,n,m)
        ie=jekc(3,n,m)
        ic=jekc(4,n,m)
        if(ie.eq.ic) then
          ivjekc(n,m)=ie
          is=ic+1
          vijekc(n,m)=-0.5d0*(fjekc(16,n,m)*(xc(is)-sx)+
     .                  0.5d0*fjekc(31,n,m)*(xc(is)-sx)**2.0d0)
          iujekc(n,m)=ic+1
          is=ie
          uijekc(n,m)=0.5d0*(fjekc(13,n,m)*(xe(is)-sx)+
     .                 0.5d0*fjekc(25,n,m)*(xe(is)-sx)**2.0d0)
        else
          ivjekc(n,m)=ie+1
          is=ic
          vijekc(n,m)=0.5d0*(fjekc(16,n,m)*(xc(is)-sx)+
     .                 0.5d0*fjekc(31,n,m)*(xc(is)-sx)**2.0d0)
          iujekc(n,m)=ic
          is=ie+1
          uijekc(n,m)=-0.5d0*(fjekc(13,n,m)*(xe(is)-sx)+
     .                  0.5d0*fjekc(25,n,m)*(xe(is)-sx)**2.0d0)
        endif
      enddo

      do n=0,njcke(m)
        sx=fjcke(3,n,m)
        ie=jcke(3,n,m)
        ic=jcke(4,n,m)
        if(ie.eq.ic) then
          iwjcke(n,m)=ie
          is=ic+1
          wijcke(n,m)=-0.5d0*(fjcke(19,n,m)*(xc(is)-sx)+
     .                  0.5d0*fjcke(37,n,m)*(xc(is)-sx)**2.0d0)
          iujcke(n,m)=ic+1
          is=ie
          uijcke(n,m)=0.5d0*(fjcke(13,n,m)*(xe(is)-sx)+
     .                 0.5d0*fjcke(25,n,m)*(xe(is)-sx)**2.0d0)
          ivjcke(n,m)=ie
          is=ic+1
          vijcke(n,m)=-0.5d0*(fjcke(16,n,m)*(xc(is)-sx)+
     .                  0.5d0*fjcke(31,n,m)*(xc(is)-sx)**2.0d0)
        else
          iwjcke(n,m)=ie+1
          is=ic
          wijcke(n,m)=0.5d0*(fjcke(19,n,m)*(xc(is)-sx)+
     .                 0.5d0*fjcke(37,n,m)*(xc(is)-sx)**2.0d0)
          iujcke(n,m)=ic
          is=ie+1
          uijcke(n,m)=-0.5d0*(fjcke(13,n,m)*(xe(is)-sx)+
     .                  0.5d0*fjcke(25,n,m)*(xe(is)-sx)**2.0d0)
          ivjcke(n,m)=ie+1
          is=ic
          vijcke(n,m)=0.5d0*(fjcke(16,n,m)*(xc(is)-sx)+
     .                 0.5d0*fjcke(31,n,m)*(xc(is)-sx)**2.0d0)
        endif
      enddo

      do n=0,njckc(m)
        sx=fjckc(3,n,m)
        ie=jckc(3,n,m)
        ic=jckc(4,n,m)
        if(ie.eq.ic) then
          iujckc(n,m)=ic+1
          is=ie
          uijckc(n,m)=0.5d0*(fjckc(13,n,m)*(xe(is)-sx)+
     .                 0.5d0*fjckc(25,n,m)*(xe(is)-sx)**2.0d0)
        else
          iujckc(n,m)=ic
          is=ie+1
          uijckc(n,m)=-0.5d0*(fjckc(13,n,m)*(xe(is)-sx)+
     .                  0.5d0*fjckc(25,n,m)*(xe(is)-sx)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end      


c-----------------------------------------------------------------------


