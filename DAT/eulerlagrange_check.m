clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat

isphr=1;
if isphr==1
  alfa1=alfa11;
else
  alfa1=alfa10;
end
alfa2=alfa20;

ns1=size(alfa1,1);
ns2=size(alfa2,1);
ms=1;

load xs0.dat
load ys0.dat
load zs0.dat

for k=1:ms
    fx(:,:,k)=xs0(1+(k-1)*ns2:k*ns2,1:ns1);
    fy(:,:,k)=ys0(1+(k-1)*ns2:k*ns2,1:ns1);
    fz(:,:,k)=zs0(1+(k-1)*ns2:k*ns2,1:ns1);
end

ks=1;
x(:,:)=fx(:,:,ks);
y(:,:)=fy(:,:,ks);
z(:,:)=fz(:,:,ks);
figure(1)
mesh(x,y,z)
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
hidden off
hold on

load ff.dat
plot3(ff(:,1),ff(:,2),ff(:,3),'ko')
hold off

load lagrange.dat;
figure(2)
mesh(alfa1,alfa2,lagrange)
xlabel('alfa1')
ylabel('alfa2')
zlabel('fs')
hold on

load euler.dat;
plot3(euler(:,1),euler(:,2),euler(:,3),'k+')
hold off

load uinterp.dat;
load vinterp.dat;
load winterp.dat;
figure(3)
mesh(alfa1,alfa2,uinterp)
xlabel('alfa1')
ylabel('alfa2')
zlabel('us')

figure(4)
mesh(alfa1,alfa2,vinterp)
xlabel('alfa1')
ylabel('alfa2')
zlabel('vs')

figure(5)
mesh(alfa1,alfa2,winterp)
xlabel('alfa1')
ylabel('alfa2')
zlabel('ws')

clear all
