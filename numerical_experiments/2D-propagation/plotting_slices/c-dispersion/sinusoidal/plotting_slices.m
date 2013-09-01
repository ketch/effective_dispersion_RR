%clear; clc; close all
figure(1); clf; hold on
set(1, 'Position', [50 50 800 300])

%% plot FV solution
FV_x_slice=load('c-dispersion_sinusoidal_xslice.txt');
FV_y_slice=load('c-dispersion_sinusoidal_yslice.txt');
xFV = FV_x_slice(:,1);
pxFV = FV_x_slice(:,2);
yFV = FV_y_slice(:,1);
pyFV = FV_y_slice(:,2);
plot(xFV,pxFV,'-b','linewidth',2)
plot(yFV,pyFV,'-r','linewidth',2)

%pause
%% plot homogenized solution
T=65; nt=T/0.5+1; %dt=0.5
%load('c-dispersion_sinusoidal.mat')
pp=squeeze(U(nt,:,:));
x=squeeze(xx(end/2,mx/2+1:mx))'-100+xFV(1); y=squeeze(yy(my/2+1:my,end/2))'-100+yFV(1);
px=pp(end/2,mx/2+1:mx); py=pp(my/2+1:my,end/2);
xi=linspace(x(1),x(end),1000); pxi=spline(x,px,xi);
yi=linspace(y(1),y(end),1000); pyi=spline(y,py,yi);
plot(xi,pxi,'--k','linewidth',2)
plot(yi,pyi,'-.k','linewidth',2)

%% formatting
xlim([60, 100])
%ylim([-0.25,0.3])
title(['t=' num2str((nt-1)*0.5)],'fontsize',20);
xlabel('x, y','fontsize',20)
leg=legend('FV x-slice', 'FV y-slice','Homogenized x-slice','Homogenized y-slice','location','northwest');
set(leg,'fontsize',20)
set(gca,'fontsize',20)

