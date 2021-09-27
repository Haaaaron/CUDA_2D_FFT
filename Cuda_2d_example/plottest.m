clear;
clc;
close all
nx=256;
ny=256;

ppp=importdata('psionek.txt');

ccc=0;
for i=1:nx
    for j=1:ny/2+1
        ccc=ccc+1;
        psik(i,j)=complex(ppp(ccc,1),ppp(ccc,2));
    end
end
psik=psik/(nx*ny);
figure; 
% pcolor((0:ny/2),(0:nx-1),sqrt(real(psik).^2+imag(psik).^2))
pcolor((0:ny/2),(0:nx-1),sqrt(real(psik.*conj(psik))))
axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

clear ppp;
ppp=importdata('psioneinv.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1),(0:nx-1),psi)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

ppp=importdata('psionekcosx.txt');

ccc=0;
for i=1:nx
    for j=1:ny/2+1
        ccc=ccc+1;
        psik(i,j)=complex(ppp(ccc,1),ppp(ccc,2));
    end
end
psik=psik/(nx*ny);
figure; 
% pcolor((0:ny/2),(0:nx-1),sqrt(real(psik).^2+imag(psik).^2))
pcolor((0:ny/2),(0:nx-1),sqrt(real(psik.*conj(psik))))
axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

clear ppp;
ppp=importdata('psioneinvcosx.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1),(0:nx-1),psi)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

ppp=importdata('psionekcosy.txt');

ccc=0;
for i=1:nx
    for j=1:ny/2+1
        ccc=ccc+1;
        psik(i,j)=complex(ppp(ccc,1),ppp(ccc,2));
    end
end
psik=psik/(nx*ny);
figure; 
% pcolor((0:ny/2),(0:nx-1),sqrt(real(psik).^2+imag(psik).^2))
pcolor((0:ny/2),(0:nx-1),sqrt(real(psik.*conj(psik))))
axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

clear ppp;
ppp=importdata('psioneinvcosy.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1),(0:nx-1),psi)
% axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar



ppp=importdata('psitwok.txt');

ccc=0;
for i=1:nx
    for j=1:ny/2+1
        ccc=ccc+1;
        psik(i,j)=complex(ppp(ccc,1),ppp(ccc,2));
    end
end
% psik=psik;
figure; 
% pcolor((0:ny/2),(0:nx-1),sqrt(real(psik).^2+imag(psik).^2))
pcolor((0:ny/2),(0:nx-1),sqrt(real(psik.*conj(psik))))
axis([0 50 0 50])
shading interp
colormap(jet)
axis square
colorbar

clear ppp;
ppp=importdata('psitwoinv.txt');

ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;
        psi(i,j)=ppp(ccc);
    end
end
figure
pcolor((0:ny-1),(0:nx-1),psi)
% axis([0 50 0 50])
% shading interp
colormap(jet)
axis square
colorbar

xxx=(1:nx)*0.8;
yyy=(1:ny)*0.8;
qt=(sqrt(3.0)/2.0)
psi(1:nx,1:ny)=0.0;
for i=1:nx
    for j=1:ny
        psi(i,j)=(cos(qt*xxx(i))*cos(qt*yyy(j)/sqrt(3.0))-0.5*cos(2*qt*yyy(j)/sqrt(3.0)));
    end
end

figure
pcolor(xxx,yyy,psi)
shading interp
colormap(jet)
axis square
colorbar
