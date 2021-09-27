clear;
clc
close all;
Angle= 0.398545
%start
nx= 2011;   ny= 3484;
Lx= 905.403526;
Ly= 1568.204909;
dx=     0.4502255227971520;
dy=     0.4501162195699461;
%stop

atone=     7.2551974569368713 ; qtone=     0.8660254037844386;
attwo=     7.2551974569368713 ;  qttwo=     0.8660254037844386;
dh_not=     9.8399999999999999 ;
deltaB_1=    -0.1500000000000000;
deltaB_2=    -0.1500000000000000;
Bx_1=     1.0000000000000000;
Bx_2=     1.0000000000000000;
tau_1=     0.8748177652797066;
tau_2=     0.8748177652797066;
v_1=     1.0000000000000000;
v_2=     1.0000000000000000;
kappa_1=     0.2081000000000000;
kappa_2=     0.2081000000000000;

dkx=2*pi/(nx*dx);
dky=2*pi/(ny*dy);

for i=1:nx
    if(i<=nx/2)
        kx(i)=(i-1)*dkx;
    else
        kx(i)=(i-1-nx)*dkx;
    end
end
for i=1:ny
    if(i<=ny/2)
        ky(i)=(i-1)*dky;
    else
        ky(i)=(i-1-ny)*dky;
    end
end
a=7.0;
for j=1:ny
    for i=1:nx
       ksq=kx(i)*kx(i)+ky(j)*ky(j);
       filtermat(i,j)=exp(-a*ksq);%ksq; %exp(-a*ksq); %ksq; %exp(-a*ksq);
    end
end



PHI=(tau_1+sqrt(tau_1^2-15.0*deltaB_1))/(15*v_1)*1.031

AAA=(9*Bx_1*PHI^2)/2;
Q=2*pi/Lx;
epps=-2.0*sin(Angle*(pi/180)) %(1.0-47.0/46.0)
H=sqrt(-2.0*(2*epps/(Q^2)+kappa_1/AAA))
% kappac=-2*epps*AAA/Q^2
% astop
 fnamepsione=dir('psione*.h5');
 fnamepsitwo=dir('psitwo*.h5');
 fnameeneloc=dir('eneloc*.h5');
 ffm='Ang00.830440FirstRun.jpg';

 nsteps=size(fnamepsitwo,1);

figure2=figure('visible', 'off');


set(gcf,'PaperPositionMode','auto')
set(figure2, 'Position', [60 100 700 600])
 ssttaa=importdata('stats.txt');
 plot(ssttaa(1:size(ssttaa,1),1),ssttaa(1:size(ssttaa,1),3),'-o');
%    print('-djpeg99','-r200',['../../../Ene',ffm])
print('-djpeg99','-r200',['Ene',ffm])
 figure1= figure('visible', 'off');
 set(gcf,'PaperPositionMode','auto')
 set(figure1, 'Position', [760 100 700 600])

for ifile=nsteps:nsteps%1:nsteps%1:1 %:nsteps % nsteps:nsteps% 1:nsteps
    fnamepsione(ifile).name
%     stop
     clf
     axes11=axes('Parent',figure1,'Position',[0.04045679012345679 0.526265060240964 0.447037037037037 0.447037037037037]);
    clear ppp;
%     ppp=importdata(fnamepsione(ifile).name);
%     ppp=importdata('psione.in');
       ppp=h5read(fnamepsione(ifile).name,'/none');
       psione=reshape(ppp,[ny,nx]);
%     stop
%     ccc=0;
%     for i=1:nx
%         for j=1:ny
%             ccc=ccc+1;klk
%             psione(i,j)=ppp(ccc);
%         end
%     end
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,psione')
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
%     ppp=importdata('psitwo.in');
    ppp=h5read(fnamepsitwo(ifile).name,'/ntwo');
    psitwo=reshape(ppp,[ny,nx]);

    ppp=h5read(fnameeneloc(ifile).name,'/ene');
    eneloc=reshape(ppp,[ny,nx]);
%     ccc=0;
%     for i=1:nx
%         for j=1:ny
%             ccc=ccc+1;
%             psitwo(i,j)=ppp(ccc);
%         end
%     end
    pcolor((0:ny-1)*dy,(0:nx-1)*dx,psitwo')
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
    m_pattern=psitwo'.*psione';
    %B = imgaussfilt(m_pattern,10,'Padding','circular');
    C=fft2(m_pattern);
    C=C.*filtermat;
    B=real(ifft2(C));

    pcolor((0:ny-1)*dy,(0:nx-1)*dx,B(1:nx-0,1:ny-0))
    shading interp
    colormap(jet)
    axis square
     axis([Ly/2-750 Ly/2+750 Lx/2-750 Lx/2+750 ])
%     axis([Ly/2-420 Ly/2 Lx/2 Lx/2+420 ])
%      axis([Ly/2-300 Ly/2+300 Lx/2-300 Lx/2+300 ])
    %axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    axis off
%     title(['Present Moire pattern','$n_2(\vec r)-n_1(\vec r)$','Interpreter','latex']);
    title(['Present ','$n_2(\vec r) \cdot n_1(\vec r)$'],'Interpreter','latex');



%     title("Present Moire pattern");

    axes22=axes('Parent',figure1,'Position',[0.5345679012345679 0.046265060240964 0.447037037037037 0.447037037037037]);

    m_pattern=eneloc';
  %  B = imgaussfilt(m_pattern,6);

   C=fft2(m_pattern);
   C=C.*filtermat;
   B=real(ifft2(C));
   %  m_pattern=psione';
   %  B = imgaussfilt(m_pattern,6);
    pcolor((40:ny-41)*dy,(40:nx-41)*dx,B(41:nx-40,41:ny-40))
    shading interp
    colormap(jet)
    axis square
     axis([Ly/2-480 Ly/2-40 Lx/2-280 Lx/2+150 ])
  %     axis([Ly/2-420 Ly/2 Lx/2 Lx/2+420 ])
  %      axis([Ly/2-300 Ly/2+300 Lx/2-300 Lx/2+300 ])
    %axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar
    axis off
  %     title(['Present Moire pattern','$n_2(\vec r)-n_1(\vec r)$','Interpreter','latex']);
    title(['Present Local energy'],'Interpreter','latex');

    shading interp
     shading flat
    colormap(jet)
    axis square
%     axis([Ly/2-280 Ly/2+280 Lx/2-280 Lx/2+280 ])
%     axis([0 Ly-dy -(Ly-Lx)/2 Lx+(Ly-Lx)/2-dx])
    colorbar

%     title([num2str(sum(sum(psione))/(nx*ny))]);
%     axis off
    fname=(['all',num2str(ifile-1,'%05.d'),'.jpg']);

%     stop
%     print('-djpeg99','-r200',fname)
%    print('-djpeg99','-r200',['../../../All', ffm])
    print('-djpeg99','-r200',['All', ffm])
    close all
    stop

end
