c-----------------------------------------------------------------------
c
      subroutine poisson2d_torus(ni,nj,d1,d2,rhs)
      include 'parameter.inc'
      include 'fft.inc'
      integer ni,nj,id,jd
      real*8 d1,d2,beta,pai,fac
      real*8 rhs(ni,nj)

      parameter(pai=
     .3.1415926535897932384626433832795028841971693993751058209749446D0)

      beta=d1/d2

      do j=1,nj
        do i=1,ni
          r1(j,i)=d1*d1*rhs(i,j)
        enddo
      enddo

      call vrfftf(nj,ni,r1,w1,nj,wsave1)

      do j=1,nj
        do i=1,ni
          r2(i,j)=r1(j,i)
        enddo
      enddo

      call vrfftf(ni,nj,r2,w2,ni,wsave2)

      id=int(ni/2)
      jd=int(nj/2)
      if((ni-id*2).eq.1) id=id+1
      if((nj-jd*2).eq.1) jd=jd+1 

      r2(1,1)=0.0d0

      do j=1,jd-1
        fac=0.5d0/(beta*beta*(dcos(dble(2*j)*pai/dble(nj))-1.0d0))
        r2(1,2*j)=fac*r2(1,2*j)
        r2(1,2*j+1)=fac*r2(1,2*j+1)
        fac=0.5d0/
     .      (beta*beta*(dcos(dble(2*j)*pai/dble(nj))-1.0d0)-2.0d0)
        r2(ni,2*j)=fac*r2(ni,2*j)
        r2(ni,2*j+1)=fac*r2(ni,2*j+1)
      enddo
      if(nj.eq.2*jd) r2(1,nj)=-0.25d0*r2(1,nj)/(beta*beta)

      do i=1,id-1
        fac=0.5/(dcos(dble(2*i)*pai/dble(ni))-1.0d0)
        r2(2*i,1)=fac*r2(2*i,1)
        r2(2*i+1,1)=fac*r2(2*i+1,1)
        fac=0.5d0/
     .      (dcos(dble(2*i)*pai/dble(ni))-1.0d0-2.0d0*beta*beta)
        r2(2*i,nj)=fac*r2(2*i,nj)
        r2(2*i+1,nj)=fac*r2(2*i+1,nj)
      enddo
      if(ni.eq.2*id) r2(ni,1)=-0.25d0*r2(ni,1)

      if(ni.eq.2*id.and.nj.eq.2*jd) then
        r2(ni,nj)=-0.25d0*r2(ni,nj)/(1.0d0+beta*beta)
      endif

      do j=1,jd-1
        do i=1,id-1
            fac=0.5d0/(dcos(dble(2*i)*pai/dble(ni))+beta*beta*
     .                 dcos(dble(2*j)*pai/dble(nj))-(1.0d0+beta*beta))
            r2(2*i,2*j)=fac*r2(2*i,2*j)
            r2(2*i,2*j+1)=fac*r2(2*i,2*j+1)
            r2(2*i+1,2*j)=fac*r2(2*i+1,2*j)
            r2(2*i+1,2*j+1)=fac*r2(2*i+1,2*j+1)
        enddo
      enddo

      call vrfftb(ni,nj,r2,w2,ni,wsave2)

      do j=1,nj
        do i=1,ni
          r1(j,i)=r2(i,j)
        enddo
      enddo

      call vrfftb(nj,ni,r1,w1,nj,wsave1)

      do j=1,nj
        do i=1,ni
          rhs(i,j)=r1(j,i)
        enddo
      enddo

      return
      end


c-----------------------------------------------------------------------
