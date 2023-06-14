c-----------------------------------------------------------------------
c
      subroutine qjc_pressure
      include 'parameter.inc'
      include 'surface.inc'
      real*8 signx,signy,signz

      DO m=1,ms

      do n=0,nicjc(m)
        signx=sign(1.0,ficjc(10,n,m))
        signy=sign(1.0,ficjc(11,n,m))
        signz=sign(1.0,ficjc(12,n,m))        
        ficjc(22,n,m)=signx*ficjc(22,n,m)
        ficjc(23,n,m)=signy*ficjc(23,n,m)
        ficjc(24,n,m)=signz*ficjc(24,n,m)
        ficjc(43,n,m)=signx*ficjc(43,n,m)
        ficjc(44,n,m)=ficjc(44,n,m)
        ficjc(45,n,m)=ficjc(45,n,m)
        ficjc(46,n,m)=signy*ficjc(46,n,m)
        ficjc(47,n,m)=ficjc(47,n,m)
        ficjc(48,n,m)=signz*ficjc(48,n,m)
      enddo

      do n=0,nickc(m)
        signx=sign(1.0,fickc(10,n,m))
        signy=sign(1.0,fickc(11,n,m))
        signz=sign(1.0,fickc(12,n,m))
        fickc(22,n,m)=signx*fickc(22,n,m)
        fickc(23,n,m)=signy*fickc(23,n,m)
        fickc(24,n,m)=signz*fickc(24,n,m)
        fickc(43,n,m)=signx*fickc(43,n,m)
        fickc(44,n,m)=fickc(44,n,m)
        fickc(45,n,m)=fickc(45,n,m)
        fickc(46,n,m)=signy*fickc(46,n,m)
        fickc(47,n,m)=fickc(47,n,m)
        fickc(48,n,m)=signz*fickc(48,n,m)
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        signy=sign(1.0,fjckc(11,n,m))
        signz=sign(1.0,fjckc(12,n,m))
        fjckc(22,n,m)=signx*fjckc(22,n,m)
        fjckc(23,n,m)=signy*fjckc(23,n,m)
        fjckc(24,n,m)=signz*fjckc(24,n,m)
        fjckc(43,n,m)=signx*fjckc(43,n,m)
        fjckc(44,n,m)=fjckc(44,n,m)
        fjckc(45,n,m)=fjckc(45,n,m)
        fjckc(46,n,m)=signy*fjckc(46,n,m)
        fjckc(47,n,m)=fjckc(47,n,m)
        fjckc(48,n,m)=signz*fjckc(48,n,m)
      enddo        

      ENDDO

      return
      end


c-----------------------------------------------------------------------
