clear;
clc;
close all
 Lx= 174.124739;  
 Ly= 175.929189; 
nx=256;
ny=256;  
 dx=     0.6801747615878317;  
 dy=     0.6872233929727672;

clear ppp;
ppp=importdata('psione00000000.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psione(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1)*dy,(0:nx-1)*dx,psione)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
colorbar
sum(sum(psione))/(nx*ny)



clear ppp;
ppp=importdata('Lsione00000000.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        Lpsione(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1)*dy,(0:nx-1)*dx,Lpsione)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
colorbar
(sum(sum(Lpsione))/(nx*ny))%/(0.75)

clear ppp;
ppp=importdata('psitwo00000000.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1)*dy,(0:nx-1)*dx,psi)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
colorbar

clear ppp;
ppp=importdata('matene00000000.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1)*dy,(0:nx-1)*dx,psi)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
colorbar
sum(sum(psi))/(nx*ny)