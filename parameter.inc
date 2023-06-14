c  parameter.inc
      implicit none

      integer i,j,k,l,m,n,m1,m2,i1,i2,i3,j1,j2,j3,k1,k2,k3
      integer nstart,nend,nstep,nx,ny,nz,ms,ns1,ns2,ij,ik,jk
      integer kip1,kip2,ks1,ks2,kks

      integer icfl,iread,iwrite,itout,iplot,ianimation,isingular,imove
      integer lw_ubc,le_ubc,ls_ubc,ln_ubc,lb_ubc,lt_ubc
      integer lw_vbc,le_vbc,ls_vbc,ln_vbc,lb_vbc,lt_vbc
      integer lw_wbc,le_wbc,ls_wbc,ln_wbc,lb_wbc,lt_wbc
      integer lw_pbc,le_pbc,ls_pbc,ln_pbc,lb_Pbc,lt_pbc

      real*8 pi,Re,Re1,cflc,cflv,t,t0,tout,dt,dt1,dt2,dtcfl,dtfix
      real*8 x0,y0,z0,xl,yl,zl,sl1,sl2,betay,betaz,dx,dy,dz
      real*8 dalfa10,dalfa11,dalfa20
      real*8 dx1,dx2,hdx,dy1,dy2,hdy,dz1,dz2,hdz

      parameter(pi=
     .3.1415926535897932384626433832795028841971693993751058209749446d0)

      parameter(Re=1.57d2,Re1=1.0d0/Re)
      parameter(x0=-4.0d0,y0=-5.0d0,z0=-5.0d0)
      parameter(xl=8.0d0,yl=8.0d0,zl=8.0d0)
      parameter(sl1=1.0d0*pi,sl2=2.0d0*pi)
      parameter(nx=240+1,ny=240+1,nz=240+1,ms=1,ns1=64,ns2=128)
      parameter(ij=2*nx*ny,ik=2*nx*nz,jk=2*ny*nz)
      parameter(kip1=2,kip2=2)
      parameter(ks1=ns1/kip1+1,ks2=ns2/kip2+1,kks=ks1*ks2)

      parameter(isingular=1,imove=1)
      parameter(lw_ubc=1,le_ubc=1,ls_ubc=1,ln_ubc=1,lb_ubc=1,lt_ubc=1)
      parameter(lw_vbc=1,le_vbc=1,ls_vbc=1,ln_vbc=1,lb_vbc=1,lt_vbc=1)
      parameter(lw_wbc=1,le_wbc=1,ls_wbc=1,ln_wbc=1,lb_wbc=1,lt_wbc=1)
      parameter(lw_pbc=2,le_pbc=2,ls_pbc=2,ln_pbc=2,lb_pbc=2,lt_pbc=2)

      common/run1/nstep,nstart,nend,
     .            icfl,iread,iwrite,itout,iplot,ianimation
      common/run2/cflc,cflv,t,t0,dt,dt1,dt2,tout,dtcfl,dtfix,
     .            dx,dy,dz,dx1,dy1,dz1,dx2,dy2,dz2,hdx,hdy,hdz,
     .            dalfa10,dalfa11,dalfa20,betay,betaz
