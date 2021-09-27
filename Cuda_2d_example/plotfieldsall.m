clear;
clc;
close all
nx=1440;
ny=2560;  
 Lx= 946.990884;  
 Ly= 1639.947475; 
 dx=     0.6576325583045440;  
 dy=     0.6406044824362107; 
 fnamepsione=dir('psione*.txt');
 fnamepsitwo=dir('psitwo*.txt');

 nsteps=size(fnamepsitwo,1);

figure1=figure;
set(gcf,'PaperPositionMode','auto')
set(figure1, 'Position', [260 600 700 600])
clear pppone
pppone=importdata(fnamepsione(1).name);
clear ppptwo
ppptwo=importdata(fnamepsitwo(1).name);
ccc=0;
for i=1:nx
    for j=1:ny
        ccc=ccc+1;       
        dpsiini(i,j)=ppptwo(ccc)-pppone(ccc);
    end
end

for ifile=nsteps:nsteps% 1:nsteps
    clf
    
    axes11=axes('Parent',figure1,'Position',[0.04045679012345679 0.526265060240964 0.447037037037037 0.447037037037037]);
    clear ppp;
    ppp=importdata(fnamepsione(ifile).name);
    
    ccc=0;
    for i=1:nx
        for j=1:ny
            ccc=ccc+1;
            psione(i,j)=ppp(ccc);
        end
    end
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,psione)
    shading interp
    colormap(jet)
    axis square
    axis([Ly/2-100 Ly/2+100 Lx/2-100 Lx/2+100 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    title('$n_1(\vec r)$','Interpreter','latex');
%     title([num2str(sum(sum(psione))/(nx*ny))]);
    axis off
%     stop
    axes12=axes('Parent',figure1,'Position',[0.5345679012345679 0.526265060240964 0.447037037037037 0.447037037037037]);
    clear ppp;
    ppp=importdata(fnamepsitwo(ifile).name);
    
    ccc=0;
    for i=1:nx
        for j=1:ny
            ccc=ccc+1;
            psitwo(i,j)=ppp(ccc);
        end
    end
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,psitwo)
    shading interp
    colormap(jet)
    axis square
    axis([Ly/2-100 Ly/2+100 Lx/2-100 Lx/2+100 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    title('$n_2(\vec r)$','Interpreter','latex');
%     title([num2str(sum(sum(psitwo))/(nx*ny))]);    
    axis off
    
    axes21=axes('Parent',figure1,'Position',[0.04045679012345679 0.046265060240964 0.447037037037037 0.447037037037037]);
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,psitwo-psione)    
    shading interp
    colormap(jet)
    axis square
    axis([Ly/2-180 Ly/2+180 Lx/2-180 Lx/2+180 ])
%      axis([Ly/2-300 Ly/2+300 Lx/2-300 Lx/2+300 ])
    %axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    axis off
%     title(['Present Moire pattern','$n_2(\vec r)-n_1(\vec r)$','Interpreter','latex']);
    title(['Present ','$n_2(\vec r)-n_1(\vec r)$'],'Interpreter','latex');
%     title("Present Moire pattern");
    
    axes22=axes('Parent',figure1,'Position',[0.5345679012345679 0.046265060240964 0.447037037037037 0.447037037037037]);   
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,dpsiini)    
    shading interp
    colormap(jet)
    axis square
    axis([Ly/2-180 Ly/2+180 Lx/2-180 Lx/2+180 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    axis off
    title(['Initial ','$n_2(\vec r)-n_1(\vec r)$'],'Interpreter','latex');
%     title(['Initial Moire pattern','$n_2(\vec r)-n_1(\vec r)$','Interpreter','latex']);
%     title("Initial Moire pattern");
    fname=(['all',num2str(ifile-1,'%05.d'),'.jpg'])
    print('-djpeg99','-r200',fname) 
    %print('-djpeg99','-r800','conftests.jpg') 
    stop
end
