clear;
clc
close all;
Angle= 0.415220    
nx= 4320;   ny= 7680; 
Lx= 2601.056589;  
 Ly= 4505.057010; 
dx=     0.6020964327507392;
dy=     0.5865959648395473; 
atone=     7.2551974569368713 ; qtone=     0.8660254037844386;
attwo=     7.2551974569368713 ;  qttwo=     0.8660254037844386; 
dh_not=    10.0700000000000003 ;  
deltaB_1=    -0.1500000000000000;  
deltaB_2=    -0.1500000000000000;  
Bx_1=     1.0000000000000000;   
Bx_2=     1.0000000000000000;  
tau_1=     0.5715000000000000;   
tau_2=     0.5715000000000000;  
v_1=     1.0000000000000000;   
v_2=     1.0000000000000000;   
kappa_1=     0.2081000000000000;   
kappa_2=     0.2081000000000000;   




PHI=(tau_1+sqrt(tau_1^2-15.0*deltaB_1))/(15*v_1)*1.031

AAA=(9*Bx_1*PHI^2)/2;
Q=2*pi/Lx;
epps=-2.0*sin(Angle*(pi/180)) %(1.0-47.0/46.0)
H=sqrt(-2.0*(2*epps/(Q^2)+kappa_1/AAA))
% kappac=-2*epps*AAA/Q^2
% astop
 fnamepsione=dir('psione*.txt');
 fnamepsitwo=dir('psitwo*.txt');
 fnamehhhone=dir('hhhone*.txt');
 fnamehhhtwo=dir('hhhtwo*.txt');

 ffm='Ang00.830440Bgg0.00001Agg1.1em1.jpg';
 
 nsteps=size(fnamepsitwo,1);
 
figure2=figure;
set(gcf,'PaperPositionMode','auto')
set(figure2, 'Position', [160 400 700 600])
 ssttaa=importdata('stats.txt');
 plot(ssttaa(2:size(ssttaa,1),1),ssttaa(2:size(ssttaa,1),3),'-o');
    print('-djpeg99','-r200',['../../../Ene',ffm]) 
 figure1=figure;
 set(gcf,'PaperPositionMode','auto')
 set(figure1, 'Position', [860 400 700 600])
% clear pppone
% pppone=importdata(fnamepsione(1).name);
% clear ppptwo
% ppptwo=importdata(fnamepsitwo(1).name);
% ccc=0;
% for i=1:nx
%     for j=1:ny
%         ccc=ccc+1;       
%         dpsiini(i,j)=ppptwo(ccc)-pppone(ccc);
%     end
% end

for ifile=nsteps:nsteps%1:nsteps%1:1 %:nsteps % nsteps:nsteps% 1:nsteps
    fnamepsione(ifile).name
     clf
     axes11=axes('Parent',figure1,'Position',[0.04045679012345679 0.526265060240964 0.447037037037037 0.447037037037037]);
    clear ppp;
%     ppp=importdata(fnamepsione(ifile).name);
    ppp=importdata('psione.in');
    
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
    axis([Ly/2-120 Ly/2+120 Lx/2-120 Lx/2+120 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    title('$n_1(\vec r)$','Interpreter','latex');
%     title([num2str(sum(sum(psione))/(nx*ny))]);
    axis off
%     stop
    axes12=axes('Parent',figure1,'Position',[0.5345679012345679 0.526265060240964 0.447037037037037 0.447037037037037]);
    clear ppp;
%     ppp=importdata(fnamepsitwo(ifile).name);
    ppp=importdata('psitwo.in');
    
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
    axis([Ly/2-120 Ly/2+120 Lx/2-120 Lx/2+120 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    title('$n_2(\vec r)$','Interpreter','latex');
%     title([num2str(sum(sum(psitwo))/(nx*ny))]);    
    axis off
    
    axes21=axes('Parent',figure1,'Position',[0.04045679012345679 0.046265060240964 0.447037037037037 0.447037037037037]);
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,-psitwo+psione)    
    shading interp
    colormap(jet)
    axis square
    axis([Ly/2-420 Ly/2 Lx/2 Lx/2+420 ])
%      axis([Ly/2-300 Ly/2+300 Lx/2-300 Lx/2+300 ])
    %axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    axis off
%     title(['Present Moire pattern','$n_2(\vec r)-n_1(\vec r)$','Interpreter','latex']);
    title(['Present ','$n_2(\vec r)-n_1(\vec r)$'],'Interpreter','latex');
%     title("Present Moire pattern");
    
    axes22=axes('Parent',figure1,'Position',[0.5345679012345679 0.046265060240964 0.447037037037037 0.447037037037037]);   
    clear hhh;
%    hhh=importdata(fnamehhhone(ifile).name);
    hhh=importdata('hhhone.in');
    
    ccc=0;
    for i=1:nx
        for j=1:ny
            ccc=ccc+1;
            hhhone(i,j)=hhh(ccc);
        end
    end
    clear hhh;
  %  hhh=importdata(fnamehhhtwo(ifile).name);
    hhh=importdata('hhhtwo.in');
    
    ccc=0;
    for i=1:nx
        for j=1:ny
            ccc=ccc+1;
            hhhtwo(i,j)=hhh(ccc);
        end
    end
%     surf((0:ny-1)*dy,(0:nx-1)*dx,hhhone); axis([Ly/2-200 Ly/2+200 Lx/2-200 Lx/2+200 ])
%     surf((0:ny-1)*dy,(0:nx-1)*dx,hhhone,hhhone)
    surf((0:ny-1)*dy,(0:nx-1)*dx,hhhone-0.5*dh_not,hhhone)
    hold on;
    surf((0:ny-1)*dy,(0:nx-1)*dx,hhhtwo+0.5*dh_not,hhhtwo)
    shading interp
     shading flat
    colormap(jet)
    axis square
%     axis([Ly/2-280 Ly/2+280 Lx/2-280 Lx/2+280 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    title('$h_{1,2}(\vec r)$','Interpreter','latex');
%     title([num2str(sum(sum(psione))/(nx*ny))]);
%     axis off
    fname=(['all',num2str(ifile-1,'%05.d'),'.jpg']);
    [max(max(hhhone)) min(min(hhhone))]
    [max(max(hhhtwo)) min(min(hhhtwo))]
%     stop
%     print('-djpeg99','-r200',fname) 
    print('-djpeg99','-r200',['../../../All', ffm]) 
%     close all
    
    figure
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,hhhone)
    shading interp
     shading flat
    colormap(jet)
    axis square
    axis([Ly/2-1300 Ly/2+1300 Lx/2-1300 Lx/2+1300 ])
    colorbar
    print('-djpeg99','-r200',['../../../hone', ffm]) 
    figure
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,hhhtwo)
    shading interp
     shading flat
    colormap(jet)
    axis square
    axis([Ly/2-1300 Ly/2+1300 Lx/2-1300 Lx/2+1300 ])
    colorbar
    print('-djpeg99','-r200',['../../../htwo',ffm]) 
    
    figure
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,-psitwo+psione)
    shading interp
     shading flat
    colormap(jet)
    axis square
    axis([Ly/2-1300 Ly/2+1300 Lx/2-1300 Lx/2+1300 ])
    colorbar
    print('-djpeg99','-r200',['../../../dpsi',ffm]) 
    clear
%     clc
    stop
end
