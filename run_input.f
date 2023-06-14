c----------------------------------------------------------------------
c
      subroutine run_input
      include 'parameter.inc'
      include 'surface.inc'

      open(unit=8,file='DAT/input.run',status='old')
      rewind 8
      read(8,*)
      read(8,*)nstep,icfl,iread,iwrite,itout,iplot,ianimation
      read(8,*)
      read(8,*)
      read(8,*)cflc,cflv,dtcfl,dtfix,tout
      close(8)

      open(unit=88,file='DAT/object.run',status='old')
      rewind 88
      read(88,*)
      read(88,*)
      read(88,*)
      read(88,*)
      read(88,*)
      do l=1,ms
        read(88,*)
        read(88,*)
        read(88,*)i01(l),ea(l),eb(l),ec(l),
     .            xsc0(l),ysc0(l),zsc0(l),
     .            phasex(l),phasey(l),phasez(l),
     .            phi0(l),theta0(l),psi0(l),
     .            phasephi(l),phasetheta(l),phasepsi(l)
      enddo
      close(88)

      open(unit=61,file='DAT/fbyfs.dat',status='unknown')
      open(unit=62,file='DAT/fbypp.dat',status='unknown')
      open(unit=63,file='DAT/xst.dat',status='unknown')
      open(unit=64,file='DAT/yst.dat',status='unknown')
      open(unit=65,file='DAT/zst.dat',status='unknown')
      open(unit=66,file='DAT/qt.dat',status='unknown')

      return
      end


c-----------------------------------------------------------------------
