clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat

ishpr=1;
ms=1;
iplot=9;
iall=1;
ilice=4;
jlice=2;
klice=2;

if ishpr==1
  alfa1=alfa11;
else
  alfa1=alfa10;
end
alfa2=alfa20;

load xc.dat
load yc.dat
load zc.dat

load xs.dat
load ys.dat
load zs.dat

xe=xc+(xc(2)-xc(1))/2;
ye=yc+(yc(2)-yc(1))/2;
ze=zc+(zc(2)-zc(1))/2;

ns1=size(alfa1,1);
ns2=size(alfa2,1);
nx=size(xc,1);
ny=size(yc,1);
nz=size(zc,1);

is=round((nx-1)/ilice);
js=round((ny-1)/jlice);
ks=round((nz-1)/klice);

for k=1:ms
  fxs(:,:,k)=xs(1+(k-1)*ns2:k*ns2,1:ns1);
  fys(:,:,k)=ys(1+(k-1)*ns2:k*ns2,1:ns1);
  fzs(:,:,k)=zs(1+(k-1)*ns2:k*ns2,1:ns1);
end

if iplot==1
  label='u';
  load u.dat
  for k=1:nz
    f(:,:,k)=u(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==2
  label='v';
  load v.dat
  for k=1:nz
    f(:,:,k)=v(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==3
  label='w';
  load w.dat
  for k=1:nz
    f(:,:,k)=w(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==4
  label='o1';
  load o1.dat
  for k=1:nz
    f(:,:,k)=o1(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==5
  label='o2';
  load o2.dat
  for k=1:nz
    f(:,:,k)=o2(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==6
  label='o3';
  load o3.dat
  for k=1:nz
    f(:,:,k)=o3(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==7
  label='p';
  load p.dat
  for k=1:nz
    f(:,:,k)=p(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==8
  label='d';
  load d.dat
  for k=1:nz
    f(:,:,k)=d(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
elseif iplot==9
  label='q';
  load q.dat
  for k=1:nz
    f(:,:,k)=q(1+(k-1)*ny:k*ny,1:nx);
  end
  fi(:,:)=f(:,is,:);
  fj(:,:)=f(js,:,:);
  fk(:,:)=f(:,:,ks);
end

figure(1)
c=ones(ns2,ns1);
for k=1:ms
  surf(fxs(:,:,k),fys(:,:,k),fzs(:,:,k),c)
  shading interp
  hold on
end
clevel=1;
piso=patch(isosurface(xc,yc,zc,f,clevel));
isonormals(xc,yc,zc,f,piso);
set(piso,'FaceColor','red','EdgeColor','none');
daspect([1 1 1]);
view(3);
camlight;
lighting phong;
xlabel('xc')
ylabel('yc')
zlabel('zc')
grid on
hold off

if iall==1

  figure(2)
  surf(zc,yc,fi)
  xlabel('zc')
  ylabel('yc')
  zlabel(label)

  figure(3)
  surf(zc,xc,fj)
  xlabel('zc')
  ylabel('xc')
  zlabel(label)

  figure(4)
  surf(xc,yc,fk)
  xlabel('xc')
  ylabel('yc')
  zlabel(label)

end

clear all

