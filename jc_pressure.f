c-----------------------------------------------------------------------
c
      subroutine jc_pressure
      include 'parameter.inc'
      include 'surface.inc'
      real*8 rp,rpp,rpm,signx,signy,signz
      real*8 duxp,duyp,duzp,dvxp,dvyp,dvzp,dwxp,dwyp,dwzp
      real*8 duxm,duym,duzm,dvxm,dvym,dvzm,dwxm,dwym,dwzm

      DO m=1,ms

      do n=0,nicjc(m)
        signx=sign(1.0,ficjc(10,n,m))
        signy=sign(1.0,ficjc(11,n,m))
        signz=sign(1.0,ficjc(12,n,m))        
        ficjc(22,n,m)=signx*ficjc(22,n,m)
        ficjc(23,n,m)=signy*ficjc(23,n,m)
        ficjc(24,n,m)=signz*ficjc(24,n,m)
        duxp=ficjc(59,n,m)
        duxm=ficjc(60,n,m)
        duyp=ficjc(61,n,m)
        duym=ficjc(62,n,m)
        duzp=ficjc(63,n,m)
        duzm=ficjc(64,n,m)
        dvxp=ficjc(65,n,m)
        dvxm=ficjc(66,n,m)
        dvyp=ficjc(67,n,m)
        dvym=ficjc(68,n,m)
        dvzp=ficjc(69,n,m)
        dvzm=ficjc(70,n,m)
        dwxp=ficjc(71,n,m)
        dwxm=ficjc(72,n,m)
        dwyp=ficjc(73,n,m)
        dwym=ficjc(74,n,m)
        dwzp=ficjc(75,n,m)
        dwzm=ficjc(76,n,m)
        rpp=duxp*dvyp-dvxp*duyp+duxp*dwzp-dwxp*duzp+dvyp*dwzp-dwyp*dvzp
        rpm=duxm*dvym-dvxm*duym+duxm*dwzm-dwxm*duzm+dvym*dwzm-dwym*dvzm
        rp=2.0d0*signz*(rpp-rpm)
        ficjc(43,n,m)=signx*(ficjc(43,n,m)+rp*ficjc(50,n,m))
        ficjc(44,n,m)=ficjc(44,n,m)+rp*ficjc(51,n,m)
        ficjc(45,n,m)=ficjc(45,n,m)+rp*ficjc(52,n,m)
        ficjc(46,n,m)=signy*(ficjc(46,n,m)+rp*ficjc(53,n,m))
        ficjc(47,n,m)=ficjc(47,n,m)+rp*ficjc(54,n,m)
        ficjc(48,n,m)=signz*(ficjc(48,n,m)+rp*ficjc(55,n,m))
      enddo

      do n=0,nickc(m)
        signx=sign(1.0,fickc(10,n,m))
        signy=sign(1.0,fickc(11,n,m))
        signz=sign(1.0,fickc(12,n,m))
        fickc(22,n,m)=signx*fickc(22,n,m)
        fickc(23,n,m)=signy*fickc(23,n,m)
        fickc(24,n,m)=signz*fickc(24,n,m)
        duxp=fickc(59,n,m)
        duxm=fickc(60,n,m)
        duyp=fickc(61,n,m)
        duym=fickc(62,n,m)
        duzp=fickc(63,n,m)
        duzm=fickc(64,n,m)
        dvxp=fickc(65,n,m)
        dvxm=fickc(66,n,m)
        dvyp=fickc(67,n,m)
        dvym=fickc(68,n,m)
        dvzp=fickc(69,n,m)
        dvzm=fickc(70,n,m)
        dwxp=fickc(71,n,m)
        dwxm=fickc(72,n,m)
        dwyp=fickc(73,n,m)
        dwym=fickc(74,n,m)
        dwzp=fickc(75,n,m)
        dwzm=fickc(76,n,m)
        rpp=duxp*dvyp-dvxp*duyp+duxp*dwzp-dwxp*duzp+dvyp*dwzp-dwyp*dvzp
        rpm=duxm*dvym-dvxm*duym+duxm*dwzm-dwxm*duzm+dvym*dwzm-dwym*dvzm
        rp=signy*2.0d0*(rpp-rpm)
        fickc(43,n,m)=signx*(fickc(43,n,m)+rp*fickc(50,n,m))
        fickc(44,n,m)=fickc(44,n,m)+rp*fickc(51,n,m)
        fickc(45,n,m)=fickc(45,n,m)+rp*fickc(52,n,m)
        fickc(46,n,m)=signy*(fickc(46,n,m)+rp*fickc(53,n,m))
        fickc(47,n,m)=fickc(47,n,m)+rp*fickc(54,n,m)
        fickc(48,n,m)=signz*(fickc(48,n,m)+rp*fickc(55,n,m))
      enddo

      do n=0,njckc(m)
        signx=sign(1.0,fjckc(10,n,m))
        signy=sign(1.0,fjckc(11,n,m))
        signz=sign(1.0,fjckc(12,n,m))
        fjckc(22,n,m)=signx*fjckc(22,n,m)
        fjckc(23,n,m)=signy*fjckc(23,n,m)
        fjckc(24,n,m)=signz*fjckc(24,n,m)
        duxp=fjckc(59,n,m)
        duxm=fjckc(60,n,m)
        duyp=fjckc(61,n,m)
        duym=fjckc(62,n,m)
        duzp=fjckc(63,n,m)
        duzm=fjckc(64,n,m)
        dvxp=fjckc(65,n,m)
        dvxm=fjckc(66,n,m)
        dvyp=fjckc(67,n,m)
        dvym=fjckc(68,n,m)
        dvzp=fjckc(69,n,m)
        dvzm=fjckc(70,n,m)
        dwxp=fjckc(71,n,m)
        dwxm=fjckc(72,n,m)
        dwyp=fjckc(73,n,m)
        dwym=fjckc(74,n,m)
        dwzp=fjckc(75,n,m)
        dwzm=fjckc(76,n,m)
        rpp=duxp*dvyp-dvxp*duyp+duxp*dwzp-dwxp*duzp+dvyp*dwzp-dwyp*dvzp
        rpm=duxm*dvym-dvxm*duym+duxm*dwzm-dwxm*duzm+dvym*dwzm-dwym*dvzm
        rp=signx*2.0d0*(rpp-rpm)
        fjckc(43,n,m)=signx*(fjckc(43,n,m)+rp*fjckc(50,n,m))
        fjckc(44,n,m)=fjckc(44,n,m)+rp*fjckc(51,n,m)
        fjckc(45,n,m)=fjckc(45,n,m)+rp*fjckc(52,n,m)
        fjckc(46,n,m)=signy*(fjckc(46,n,m)+rp*fjckc(53,n,m))
        fjckc(47,n,m)=fjckc(47,n,m)+rp*fjckc(54,n,m)
        fjckc(48,n,m)=signz*(fjckc(48,n,m)+rp*fjckc(55,n,m))
      enddo        

      ENDDO

      return
      end


c-----------------------------------------------------------------------
