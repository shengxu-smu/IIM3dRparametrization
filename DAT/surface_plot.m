clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat

ms=1;
ks=1;

plot=1;
ishpr=1;
if ishpr==1
  alfa1=alfa11;
else
  alfa1=alfa10;
end
alfa2=alfa20;

ns1=size(alfa1,1);
ns2=size(alfa2,1);


if plot==0
  load use.dat
  load vse.dat
  load wse.dat
elseif plot==1
  load us.dat
  load vs.dat
  load ws.dat
elseif plot==-1
  load use.dat
  load vse.dat
  load wse.dat
  load us.dat
  load vs.dat
  load ws.dat
end
load fs.dat

for k=1:ms
  if plot==0
    fu(:,:,k)=use(1+(k-1)*ns2:k*ns2,1:ns1);
    fv(:,:,k)=vse(1+(k-1)*ns2:k*ns2,1:ns1);
    fw(:,:,k)=wse(1+(k-1)*ns2:k*ns2,1:ns1);
  elseif plot==1
    fu(:,:,k)=us(1+(k-1)*ns2:k*ns2,1:ns1);
    fv(:,:,k)=vs(1+(k-1)*ns2:k*ns2,1:ns1);
    fw(:,:,k)=ws(1+(k-1)*ns2:k*ns2,1:ns1);
  elseif plot==-1
    fu(:,:,k)=us(1+(k-1)*ns2:k*ns2,1:ns1)-use(1+(k-1)*ns2:k*ns2,1:ns1);
    fv(:,:,k)=vs(1+(k-1)*ns2:k*ns2,1:ns1)-vse(1+(k-1)*ns2:k*ns2,1:ns1);
    fw(:,:,k)=ws(1+(k-1)*ns2:k*ns2,1:ns1)-wse(1+(k-1)*ns2:k*ns2,1:ns1);
  end
  ff(:,:,k)=fs(1+(k-1)*ns2:k*ns2,1:ns1);
end

u(:,:)=fu(:,:,ks);
v(:,:)=fv(:,:,ks);
w(:,:)=fw(:,:,ks);
f(:,:)=ff(:,:,ks);

figure(1)
mesh(alfa1,alfa2,u)
xlabel('alfa1')
ylabel('alfa2')
zlabel('us')

figure(2)
mesh(alfa1,alfa2,v)
xlabel('alfa1')
ylabel('alfa2')
zlabel('vs')

figure(3)
mesh(alfa1,alfa2,w)
xlabel('alfa1')
ylabel('alfa2')
zlabel('ws')

figure(4)
mesh(alfa1,alfa2,f)
xlabel('alfa1')
ylabel('alfa2')
zlabel('fn')

