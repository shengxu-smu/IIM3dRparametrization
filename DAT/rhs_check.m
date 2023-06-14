clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat
alfa1=alfa11;
alfa2=alfa20;

ns1=size(alfa1,1);
ns2=size(alfa2,1);

load tmp1.dat
load tmp2.dat
load tmp3.dat
load tmp4.dat
load tmp5.dat
load tmp6.dat
load tmp7.dat
load tmp8.dat

figure(1)
mesh(alfa1,alfa2,tmp1)

figure(2)
mesh(alfa1,alfa2,tmp2)

figure(3)
mesh(alfa1,alfa2,tmp3)

figure(4)
mesh(alfa1,alfa2,tmp4)

figure(5)
mesh(alfa1,alfa2,tmp5-tmp1)

figure(6)
mesh(alfa1,alfa2,tmp6-tmp2)

figure(7)
mesh(alfa1,alfa2,tmp7-tmp3)

figure(8)
mesh(alfa1,alfa2,tmp8-tmp4)

clear all
