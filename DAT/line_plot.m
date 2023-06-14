clear all

load xc.dat
load yc.dat
load zc.dat

load xu.dat
load yv.dat
load zw.dat
load xp.dat
load yp.dat
load zp.dat

load fbyfs.dat
cxcycz=fbyfs;
ks=1;
f=[cxcycz(:,1),cxcycz(:,2+(ks-1)*3:1+3*ks)];

figure(1)
plot(xc,xu,'bo-',yc,yv,'ro-',zc,zw,'go-')
legend('u(x)','v(y)','w(z)')

figure(2)
plot(xc,xp,'bo-',yc,yp,'ro-',zc,zp,'go-')
legend('p(x)','p(y)','p(z)')

tp=2*pi/0.8;
fac=6/(5+pi/4);
figure(3)
plot(f(:,1)/tp,fac*f(:,2),'b-',f(:,1)/tp,fac*f(:,3),'r-',f(:,1)/tp,fac*f(:,4),'g-')
legend('cx','cy','cz')

