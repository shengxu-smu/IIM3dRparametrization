c-----------------------------------------------------------------------
c
      subroutine correction_strain
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'
      integer ie,je,ke,ic,jc,kc,ie1,je1,ke1
      real*8 sx,sy,sz

      DO m=1,ms

      do n=0,nicjc(m)
        sz=ficjc(3,n,m)
        ke=icjc(3,n,m)
        ke1=ke+1
        kc=icjc(4,n,m)
        if(ke.eq.kc) then
          udz(n,m)=-dz1*(ficjc(15,n,m)*(ze(ke)-sz)+
     .             0.5d0*ficjc(30,n,m)*(ze(ke)-sz)**2.0d0)
          vdz(n,m)=-dz1*(ficjc(18,n,m)*(ze(ke)-sz)+
     .             0.5d0*ficjc(36,n,m)*(ze(ke)-sz)**2.0d0)
          wdz(n,m)=-dz1*(ficjc(21,n,m)*(ze(ke)-sz)+
     .             0.5d0*ficjc(42,n,m)*(ze(ke)-sz)**2.0d0)
        else
          udz(n,m)=-dz1*(ficjc(15,n,m)*(ze(ke1)-sz)+
     .             0.5d0*ficjc(30,n,m)*(ze(ke1)-sz)**2.0d0)
          vdz(n,m)=-dz1*(ficjc(18,n,m)*(ze(ke1)-sz)+
     .             0.5d0*ficjc(36,n,m)*(ze(ke1)-sz)**2.0d0)
          wdz(n,m)=-dz1*(ficjc(21,n,m)*(ze(ke1)-sz)+
     .             0.5d0*ficjc(42,n,m)*(ze(ke1)-sz)**2.0d0)
        endif
      enddo

      do n=0,nickc(m)
        sy=fickc(3,n,m)
        je=ickc(3,n,m)
        je1=je+1
        jc=ickc(4,n,m)
        if(je.eq.jc) then
          udy(n,m)=-dy1*(fickc(14,n,m)*(ye(je)-sy)+
     .             0.5d0*fickc(28,n,m)*(ye(je)-sy)**2.0d0)
          vdy(n,m)=-dy1*(fickc(17,n,m)*(ye(je)-sy)+
     .             0.5d0*fickc(34,n,m)*(ye(je)-sy)**2.0d0)
          wdy(n,m)=-dy1*(fickc(20,n,m)*(ye(je)-sy)+
     .             0.5d0*fickc(40,n,m)*(ye(je)-sy)**2.0d0)
        else
          udy(n,m)=-dy1*(fickc(14,n,m)*(ye(je1)-sy)+
     .             0.5d0*fickc(28,n,m)*(ye(je1)-sy)**2.0d0)
          vdy(n,m)=-dy1*(fickc(17,n,m)*(ye(je1)-sy)+
     .             0.5d0*fickc(34,n,m)*(ye(je1)-sy)**2.0d0)
          wdy(n,m)=-dy1*(fickc(20,n,m)*(ye(je1)-sy)+
     .             0.5d0*fickc(40,n,m)*(ye(je1)-sy)**2.0d0)
        endif
      enddo

      do n=0,njckc(m)
        sx=fjckc(3,n,m)
        ie=jckc(3,n,m)
        ie1=ie+1
        ic=jckc(4,n,m)
        if(ie.eq.ic) then
          udx(n,m)=-dx1*(fjckc(13,n,m)*(xe(ie)-sx)+
     .             0.5d0*fjckc(25,n,m)*(xe(ie)-sx)**2.0d0)
          vdx(n,m)=-dx1*(fjckc(16,n,m)*(xe(ie)-sx)+
     .             0.5d0*fjckc(31,n,m)*(xe(ie)-sx)**2.0d0)
          wdx(n,m)=-dx1*(fjckc(19,n,m)*(xe(ie)-sx)+
     .             0.5d0*fjckc(37,n,m)*(xe(ie)-sx)**2.0d0)
        else
          udx(n,m)=-dx1*(fjckc(13,n,m)*(xe(ie1)-sx)+
     .             0.5d0*fjckc(25,n,m)*(xe(ie1)-sx)**2.0d0)
          vdx(n,m)=-dx1*(fjckc(16,n,m)*(xe(ie1)-sx)+
     .             0.5d0*fjckc(31,n,m)*(xe(ie1)-sx)**2.0d0)
          wdx(n,m)=-dx1*(fjckc(19,n,m)*(xe(ie1)-sx)+
     .             0.5d0*fjckc(37,n,m)*(xe(ie1)-sx)**2.0d0)
        endif
      enddo

      ENDDO

      return
      end      


c-----------------------------------------------------------------------


