clear all

load xc.dat
load yc.dat
load zc.dat

nx=size(xc,1);
ny=size(yc,1);
nz=size(zc,1);

lice=2;
is=round((nx-1)/lice);
js=round((ny-1)/lice);
ks=round((nz-1)/lice);

iplot=4;
if iplot==1
  label='ue';
  load ue.dat
  for k=1:nz
    f(:,:,k)=ue(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==2
  label='ve';
  load ve.dat
  for k=1:nz
    f(:,:,k)=ve(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==3
  label='we';
  load we.dat
  for k=1:nz
    f(:,:,k)=we(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==4
  label='pe';
  load pe.dat
  for k=1:nz
    f(:,:,k)=pe(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
end

figure(1)
surf(zc,yc,fi)
xlabel('zc')
ylabel('yc')
zlabel(label)

figure(2)
surf(zc,xc,fj)
xlabel('zc')
ylabel('xc')
zlabel(label)

figure(3)
surf(xc,yc,fk)
xlabel('xc')
ylabel('yc')
zlabel(label)

clear all
