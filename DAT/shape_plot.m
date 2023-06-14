clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat
alfa1=alfa10;
alfa2=alfa20;

ns1=size(alfa1,1);
ns2=size(alfa2,1);
ms=1;
initial=1;

if initial==1
  load xs0.dat
  load ys0.dat
  load zs0.dat
else
  load xs.dat
  load ys.dat
  load zs.dat
end

if initial==1
  for k=1:ms
    fx(:,:,k)=xs0(1+(k-1)*ns2:k*ns2,1:ns1);
    fy(:,:,k)=ys0(1+(k-1)*ns2:k*ns2,1:ns1);
    fz(:,:,k)=zs0(1+(k-1)*ns2:k*ns2,1:ns1);
  end
else
  for k=1:ms
    fx(:,:,k)=xs(1+(k-1)*ns2:k*ns2,1:ns1);
    fy(:,:,k)=ys(1+(k-1)*ns2:k*ns2,1:ns1);
    fz(:,:,k)=zs(1+(k-1)*ns2:k*ns2,1:ns1);
  end
end

ks=1;
x(:,:)=fx(:,:,ks);
y(:,:)=fy(:,:,ks);
z(:,:)=fz(:,:,ks);
c=ones(ns2,ns1);

figure(1)
mesh(x,y,z,c)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
hold on

if ms==2
  ks=2;
  x(:,:)=fx(:,:,ks);
  y(:,:)=fy(:,:,ks);
  z(:,:)=fz(:,:,ks);
  mesh(x,y,z)
  axis equal
end
hold off

figure(2)
surf(alfa11,alfa20,x)

figure(3)
surf(alfa11,alfa20,y)

figure(4)
surf(alfa11,alfa20,z)

clear all



