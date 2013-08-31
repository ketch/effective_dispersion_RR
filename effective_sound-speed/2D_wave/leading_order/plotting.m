clear; clc; 


load('sound-speed_2D-wave.mat')

% plotting surface plot at t=0
T=0; nt=T/0.5+1; %dt=0.5
pp=squeeze(U(nt,:,:));
figure(1); clf; hold on
set(1, 'Position', [50 50 700 300])
surf(xx,yy,pp,'EdgeColor','none');
colorbar; view(2)
title(['t=' num2str(0.5*(nt-1))],'fontsize',20);
axis equal
set(gca,'fontsize',20)

% plotting surface plot at t=5
T=5; nt=T/0.5+1; %dt=0.5
pp=squeeze(U(nt,:,:));
figure(2); clf; hold on
set(2, 'Position', [50 50 700 300])
surf(xx,yy,pp,'EdgeColor','none');
colorbar; view(2)
title(['t=' num2str(0.5*(nt-1))],'fontsize',20);
axis equal
set(gca,'fontsize',20)
