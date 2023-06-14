c field.inc
      integer ind(0:2*nx-3,0:2*ny-3,0:2*nz-3)

      real*8 x(-2:2*nx),y(-2:2*ny),z(-2:2*nz)
      real*8 xe(0:nx),ye(0:ny),ze(0:nz)
      real*8 xc(0:nx+1),yc(0:ny+1),zc(0:nz+1)

      real*8 pw(ny,nz),pe(ny,nz)
      real*8 ps(nx,nz),pn(nx,nz)
      real*8 pb(nx,ny),pt(nx,ny)

      real*8 u(0:nx,0:ny+1,0:nz+1)
      real*8 v(0:nx+1,0:ny,0:nz+1)
      real*8 w(0:nx+1,0:ny+1,0:nz)
      real*8 p(0:nx+1,0:ny+1,0:nz+1)      
      real*8 d(0:nx+1,0:ny+1,0:nz+1)
      real*8 o(0:nx+1,0:ny+1,0:nz+1)

      real*8 un(0:nx,0:ny+1,0:nz+1)
      real*8 vn(0:nx+1,0:ny,0:nz+1)
      real*8 wn(0:nx+1,0:ny+1,0:nz)
      real*8 dn(0:nx+1,0:ny+1,0:nz+1)

      real*8 um(0:nx,0:ny+1,0:nz+1)
      real*8 vm(0:nx+1,0:ny,0:nz+1)
      real*8 wm(0:nx+1,0:ny+1,0:nz)

      real*8 fcee(0:nx+1,0:ny+1,0:nz+1)
      real*8 fcec(0:nx+1,0:ny+1,0:nz+1)
      real*8 fece(0:nx+1,0:ny+1,0:nz+1)
      real*8 fecc(0:nx+1,0:ny+1,0:nz+1)
      real*8 fcce(0:nx+1,0:ny+1,0:nz+1)
      real*8 fccc(0:nx+1,0:ny+1,0:nz+1)
      real*8 feec(0:nx+1,0:ny+1,0:nz+1)

      real*8 data1(0:nx+1,0:ny+1,0:nz+1)
      real*8 data2(0:nx+1,0:ny+1,0:nz+1)
      real*8 data3(0:nx+1,0:ny+1,0:nz+1)
      real*8 data4(0:nx+1,0:ny+1,0:nz+1)
      real*8 data5(0:nx+1,0:ny+1,0:nz+1)
      real*8 data6(0:nx+1,0:ny+1,0:nz+1)
      real*8 data7(0:nx+1,0:ny+1,0:nz+1)
      real*8 data8(0:nx+1,0:ny+1,0:nz+1)
      real*8 data9(0:nx+1,0:ny+1,0:nz+1)

      common/index/ind
      common/grid/x,y,z,xe,ye,ze,xc,yc,zc
      common/bc/pw,pe,ps,pn,pb,pt
      common/field/u,v,w,p,d,o,un,vn,wn,dn,um,vm,wm
      common/extra/fcee,fcec,fece,fecc,fcce,fccc,feec
      common/swap/data1,data2,data3,data4,
     .            data5,data6,data7,data8,data9
