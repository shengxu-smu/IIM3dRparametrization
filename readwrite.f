c-----------------------------------------------------------------------
c
      subroutine data_read
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      open(unit=11,file='DAT/old.out',form='unformatted',status='old')
      rewind 11
      read(11)nstart,t,t0
      read(11)xsc,ysc,zsc,phi,theta,psi
      read(11)xs,ys,zs,us,vs,ws
      read(11)x,y,z,xe,ye,ze,xc,yc,zc
      read(11)u,v,w,p
      close(11)

      return
      end


c-----------------------------------------------------------------------
c
      subroutine data_write
      include 'parameter.inc'
      include 'surface.inc'
      include 'field.inc'

      open(unit=12,file='DAT/new.out',form='unformatted',
     .     status='unknown')
      rewind 12
      write(12)nend,t,t0
      write(12)xsc,ysc,zsc,phi,theta,psi
      write(12)xs,ys,zs,us,vs,ws
      write(12)x,y,z,xe,ye,ze,xc,yc,zc
      write(12)u,v,w,p
      close(12)

      return
      end


c-----------------------------------------------------------------------
