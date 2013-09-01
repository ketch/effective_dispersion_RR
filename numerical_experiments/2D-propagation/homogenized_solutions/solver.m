function []=solver()
clc; clear all; clf

%compute parameters given: K1, p1, K2, p2

%%parameters_layered(4,4,1,1); name_solution='z-dispersion.mat';
%%parameters_layered(5/8,8/5,5/2,2/5); name_solution='c-dispersion_layered.mat';

%%parameters_sinusoidal(5/8,8/5,5/2,2/5); name_solution='c-dispersion_sinusoidal.mat';

%%parameters_layered(16,1,1,1); name_solution='iso-high-dispersion.mat';
%%parameters_layered(16,1,1,1); name_solution='iso-low-dispersion.mat'; %NOTE: change variance to 5

%%parameters_layered(1,1,1.5,2.5); name_solution='cz-dispersion.mat';

%VARIANCE
varx=2; vary=2;

dt=0.01;    %time step
tf=100;      %final time
td=0.5;     %time interval to display 

save_solution = 1; %flag to save solution

% physical domain
x_lower=0; x_upper=200;
y_lower=0; y_upper=200;

mx=2^8; %Number of Fourier modes
my=2^8;
Lx=x_upper-x_lower;
kx = (2*pi/Lx)*[0:(mx/2-1) (-mx/2):-1]; % Wavenumber vector in x
Ly=(y_upper-y_lower);
ky = (2*pi/Ly)*[0:(my/2-1) (-my/2):-1]; % Wavenumber vector in x

%discretized domain
dx = (x_upper-x_lower)/mx;
dy = (y_upper-y_lower)/my;
x = (0:(mx-1))*dx;
y = (0:(my-1))*dy;
[xx, yy] = meshgrid(x,y);

nit=floor(tf/dt); %number of iterations

%initial conditions
A=5;
x0=(x_upper-x_lower)/2; y0=(y_upper-y_lower)/2;
s=A*exp(-(xx-x0).^2/(2*varx)-(yy-y0).^2/(2*vary)); %IC het in x and y
%s=A*exp(-(xx-x0).^2/(2*varx)); %IC hom in y
%s=A*exp(-(yy-y0).^2/(2*vary)); %IC hom in x

u(1,:,:)=xx.*0; %u
u(2,:,:)=xx.*0; %v
u(3,:,:)=s;   %sig

ss=squeeze(u(3,:,:));
surf(xx,yy,ss,'EdgeColor','none');
colorbar; view(2)
axis equal
pause(2)

U(1,:,:)=ss;

index=1;

for i=1:nit
    disp('*********************************')
    disp('*********************************')
    disp(['Time step ' num2str(i) '. Time t=' num2str(i*dt)])
        % Four stages Runge-Kutta
        D1u=dt.*ps_discretization(u,kx,ky);    
        D2u=dt.*ps_discretization(u+0.5*D1u,kx,ky);
        D3u=dt.*ps_discretization(u+0.5*D2u,kx,ky);
        D4u=dt.*ps_discretization(u+D3u,kx,ky);
        u = u + (D1u+2*D2u+2*D3u+D4u)/6;
    if((i*dt-index*td)>=0)
        ss=squeeze(u(3,:,:));
        surf(xx,yy,ss,'EdgeColor','none');
        colorbar; view(2)
        title(['t=' num2str(dt*i)]);
        U(index+1,:,:)=ss;
        index=index+1;
        axis equal
        pause(0.1)

    end
end
if save_solution == 1
    save(name_solution)
end
