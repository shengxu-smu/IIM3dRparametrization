c-----------------------------------------------------------------------
c
      subroutine surface_plot
      include 'parameter.inc'
      include 'surface.inc'

      write(*,*)
      write(*,*)'odd numbers:'
      write(*,*)'niejc = ',(niejc(l),l=1,ms)
      write(*,*)'nicje = ',(nicje(l),l=1,ms)
      write(*,*)'nicjc = ',(nicjc(l),l=1,ms)
      write(*,*)'niekc = ',(niekc(l),l=1,ms)
      write(*,*)'nicke = ',(nicke(l),l=1,ms)
      write(*,*)'nickc = ',(nickc(l),l=1,ms)
      write(*,*)'njekc = ',(njekc(l),l=1,ms)
      write(*,*)'njcke = ',(njcke(l),l=1,ms)
      write(*,*)'njckc = ',(njckc(l),l=1,ms)

      open(unit=70,file='DAT/fs.dat',status='unknown')
      open(unit=71,file='DAT/xs.dat',status='unknown')
      open(unit=72,file='DAT/ys.dat',status='unknown')
      open(unit=73,file='DAT/zs.dat',status='unknown')
      open(unit=74,file='DAT/us.dat',status='unknown')
      open(unit=75,file='DAT/vs.dat',status='unknown')
      open(unit=76,file='DAT/ws.dat',status='unknown')

      DO m=1,ms
      do m2=0,ns2
        write(70,200)(fr(m1,m2,4,m),m1=0,ns1)
        write(71,200)(xs(m1,m2,m),m1=0,ns1)
        write(72,200)(ys(m1,m2,m),m1=0,ns1)
        write(73,200)(zs(m1,m2,m),m1=0,ns1)
        write(74,200)(us(m1,m2,m),m1=0,ns1)
        write(75,200)(vs(m1,m2,m),m1=0,ns1)
        write(76,200)(ws(m1,m2,m),m1=0,ns1)
      enddo
      ENDDO

      close(70)
      close(71)
      close(72)
      close(73)
      close(74)
      close(75)
      close(76)

200   format(1x,1000e16.6e4)

      return
      end


c-----------------------------------------------------------------------
