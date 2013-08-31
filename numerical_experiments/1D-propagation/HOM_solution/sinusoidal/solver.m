function solver
clc; clear

t=970;
xshift=1000;

%% non-zoomed plot
[x,s0]=stress(t,0,0,0); [x,s2]=stress(t,1,0,0);
[x,s4]=stress(t,1,1,0); [x,s6]=stress(t,1,1,1);
x=x-200+xshift;
figure(1); clf; hold on
set(1, 'Position', [50 50 800 300])
set(gca,'FontSize',20)
plot(x,s0,'Color',[0.15,0.65,0.15]); 
plot(x,s2,'r')
plot(x,s4,'Color',[1,0.4,0.6]); 
plot(x,s6,'b')
title(['Sinusoidal domain at t=' num2str(t)],'fontsize',20)
xlim([50+xshift 175+xshift])
% plot finite volume solution
FV=textread('stress.txt');
xFV=FV(:,1); xFV=linspace(xshift,200+xshift,length(xFV))+xFV(1);
sFV=FV(:,2);
plot(xFV,sFV,'--k')
leg=legend('leading order','2nd correction','4th correction','6th correction','FV solution','location','northwest');
set(leg,'FontSize',20)

%% zoomed plot
[x,s0]=stress(t,0,0,0); [x,s2]=stress(t,1,0,0);
[x,s4]=stress(t,1,1,0); [x,s6]=stress(t,1,1,1);
x=x-200+xshift;
figure(2); clf; hold on
set(2, 'Position', [50 50 800 300])
set(gca,'FontSize',20)
plot(x,s0,'Color',[0.15,0.65,0.15]); 
plot(x,s2,'r')
plot(x,s4,'Color',[1,0.4,0.6]); 
plot(x,s6,'b')
title(['Sinusoidal domain at t=' num2str(t)],'fontsize',20)
axis([79+xshift 81+xshift -0.075 0.075])
% plot finite volume solution
FV=textread('stress.txt');
xFV=FV(:,1); xFV=linspace(xshift,200+xshift,length(xFV))+xFV(1); 
sFV=FV(:,2);
plot(xFV,sFV,'--k')
%
FV=textread('stress_normal.txt');
xFV=FV(:,1); xFV=linspace(xshift,200+xshift,length(xFV))+xFV(1);
sFV=FV(:,2);
plot(xFV,sFV,'--c')
leg=legend('leading order','2nd correction','4th correction','6th correction',...
    'FV finer solution','FV coarser solution','location','northwest');
set(leg,'FontSize',20)


function [x,stress] = stress(t,delta2,delta4,delta6)
%coefficients
a2 = delta2*  -0.013208894074369;
a4 = delta4*  -1.817206413939849e-04;
a6 = delta6*  6.071198284057355e-06;
b2  = delta2*   -0.011033010870322;
b4  = delta4*   -2.347404625178936e-05;
b6  = delta6*   6.906004360487817e-06;

Kh=1.25; ph=0.64;
ceff=sqrt(Kh/ph);

% physical domain
x_lower=0; x_upper=400;
mx=2^10; %Number of Fourier modes
Lx=x_upper-x_lower;
kx = (2*pi/Lx)*[0:(mx/2-1) (-mx/2):-1]; % Wavenumber vector in x

%discretized domain
dx = (x_upper-x_lower)/mx;
x = (0:(mx-1))*dx;

%initial conditions
A=5; %amp is half of FV since another half travels to the left
x0=(x_upper-x_lower)/2;
varx=5;
s=A*exp(-(x-x0).^2/(2*varx)); %IC hom in y

% velocity as a function of wave number
delta=1;
c=@(k)ceff.*sqrt(1+delta^2*(a2+b2)*k.^2 ...
        -delta^4*(-a2*b2+b4+a4)*k.^4 ...
        +delta^6*(-a2*b4-a4*b2+b6+a6)*k.^6 ...
        +delta^8*(a2*b6+a4*b4+a6*b2)*k.^8 ...
        -delta^10*(a4*b6+a6*b4)*k.^10 ...
        +delta^12*(a6*b6)*k.^12 );

% solve diferential equation
shat=fft(s);
shift=t*c(kx);
sshifted=real(ifft2(exp(-1i*kx.*shift).*shat));
stress=sshifted;

end

end