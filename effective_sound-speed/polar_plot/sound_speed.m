clear all; clc

KA=5/8; rhoA=8/5; 
KB=5/2; rhoB=2/5; 

Kh=2*KA*KB/(KA+KB);
rhoh=2*rhoA*rhoB/(rhoA+rhoB);
rhom=(rhoA+rhoB)/2;

Kh=1;
rhoh=1;
rhom=[1 2 4 8];
color=['b' 'r' 'c' 'k'];
%close all
figure(1)
clf; %hold on
for i=1:length(rhom)
    th=linspace(0,2*pi,100);
    ceff=sqrt(Kh/rhoh/rhom(i))*sqrt(rhom(i)*cos(th).^2+rhoh*sin(th).^2);
    figure(gcf)
    h1=polar(th,ceff);
    set(h1,'color',color(i),'linewidth',2)
    hold on
end
hold off
axis equal
