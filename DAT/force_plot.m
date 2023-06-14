load cdcl.dat;
load cxcycz.dat;
[m1,n1]=size(cdcl);
[m2,n2]=size(cxcycz);
tp=2*pi/0.8;
fac=4/(3+pi/4);
k1=10;
k2=10;

figure(1)
plot(cdcl(1:k1:m1,1)/tp,cdcl(1:k1:m1,2),'r-',cxcycz(1:k2:m2,1)/tp,fac*cxcycz(1:k2:m2,3),'g-')
xlabel('t')
ylabel('cd')
legend('2d','3d')

figure(2)
plot(cdcl(1:k1:m1,1)/tp,cdcl(1:k1:m1,3),'r-',cxcycz(1:k2:m2,1)/tp,fac*cxcycz(1:k2:m2,4),'g-')
xlabel('t')
ylabel('cl')
legend('2d','3d')

