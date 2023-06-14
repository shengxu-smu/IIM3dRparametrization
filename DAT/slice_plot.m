clear all

load xc.dat
load yc.dat
load zc.dat

load uslc.dat
load vslc.dat
load wslc.dat
load pslc.dat
load o1slc.dat
load o2slc.dat
load o3slc.dat

lice=1;
intvl=1;
iquiver=0;

if lice==1
  label='o1';
  m=size(yc,1);
  n=size(zc,1);
elseif lice==2
  label='o2';
  m=size(xc,1);
  n=size(zc,1);
elseif lice==3
  label='o3';
  m=size(xc,1);
  n=size(yc,1);
end

figure(1)
if lice==1  
  pcolor(yc,zc,o1slc);
  xlabel('y')
  ylabel('z')
elseif lice==2
  pcolor(xc,zc,o2slc);
  xlabel('x')
  ylabel('z')
elseif lice==3
  pcolor(xc,yc,o3slc);
  xlabel('x')
  ylabel('y')
end
title(label)
shading interp;
caxis([-5 5]);
axis equal;
axis([-1 3 -1 3])
hold on
if iquiver==1
  if lice==1
    quiver(yc(1:intvl:m),zc(1:intvl:n),...
           vslc(1:intvl:n,1:intvl:m),wslc(1:intvl:n,1:intvl:m))
  elseif lice==2
    quiver(xc(1:intvl:m),zc(1:intvl:n),...
           uslc(1:intvl:n,1:intvl:m),wslc(1:intvl:n,1:intvl:m))
  elseif lice==1
    quiver(xc(1:intvl:m),yc(1:intvl:n),...
           uslc(1:intvl:n,1:intvl:m),vslc(1:intvl:n,1:intvl:m))  
  end
end
hold off

figure(2)
label='p';
if lice==1
  pcolor(yc,zc,pslc);
  xlabel('y')
  ylabel('z')
elseif lice==2
  pcolor(xc,zc,pslc);
  xlabel('x')
  ylabel('z')
elseif lice==3
  pcolor(xc,yc,pslc);
  xlabel('x')
  ylabel('y')
end
title(label)
shading interp;
caxis([-1 1]);
axis equal;
axis([-1 3 -1 3])

clear all
