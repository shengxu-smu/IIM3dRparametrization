clear all

load alfa10.dat
load alfa11.dat
load alfa20.dat
ishpr=1;
ms=1;
if ishpr==1
  alfa1=alfa11;
else
  alfa1=alfa10;
end
alfa2=alfa20;

load xc.dat
load yc.dat
load zc.dat
ns1=size(alfa1,1);
ns2=size(alfa2,1);
nx=size(xc,1);
ny=size(yc,1);
nz=size(zc,1);

% Control no. of frames
numframes = 20;

% Control how fast final movie will play; total time
% of movie in seconds will be (numframes)/(num_frames_per_second)
num_frames_per_second = 2;
dur = numframes/num_frames_per_second;
disp('                                              ');
disp('                                              ');
disp(['     duration of movie will be  ',num2str(dur), '  secs']);

% Create a new Audio Video Interleaved (AVI) file to be called
% something.avi using the avifile function
aviobj = avifile('flapper3D.avi','fps',num_frames_per_second);

% Create plots of 3_D and add the 'frames' to the avi movie
% using the addframe function

xstid=fopen('xst.dat', 'rt');
ystid=fopen('yst.dat', 'rt');
zstid=fopen('zst.dat', 'rt');
qtid=fopen('qt.dat', 'rt');
for kk = 1:numframes
  c=ones(ns1,ns2);
  for k=1:ms
    fxs=fscanf(xstid,'%e',[ns1,ns2]);
    fys=fscanf(ystid,'%e',[ns1,ns2]);
    fzs=fscanf(zstid,'%e',[ns1,ns2]);
    surf(fxs,fys,fzs,c)
    shading interp
    axis equal
    axis([-4 4 -2 4 -2 4])
    hold on
  end

  q=fscanf(qtid,'%e',[nx,ny*nz]);
  qq=q';
  for k=1:nz
    f(:,:,k)=qq(1+(k-1)*ny:k*ny,1:nx);
  end
  clevel=1;
  piso=patch(isosurface(xc,yc,zc,f,clevel));
  isonormals(xc,yc,zc,f,piso);
  set(piso,'FaceColor','red','EdgeColor','none');
  daspect([1 1 1]);
  view(3);
  camlight;
  lighting phong;
  hold off

  frame = getframe;
  aviobj = addframe(aviobj,frame);
end
fclose(xstid);
fclose(ystid);
fclose(zstid);
fclose(qtid);

% Tells matlab that there are no more frames to be added
% and to close the 2Dsinwave.avi file
aviobj = close(aviobj);

clear all
