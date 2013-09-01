function parameters_sinusoidal(K1,p1,K2,p2)
clc
addpath '/Users/mquezada/Desktop/dispersion-relation/chebfun'
global alpha1 alpha2 alpha3 alpha4 alpha5 beta1 beta2 beta3 beta4 beta5 gamma1 gamma2 gamma3 gamma4 gamma5 Kh ph pm

recompute_coefficients=1;

if recompute_coefficients==0
    pm=0.8;
    ph=0.64;
    Kh=1.25;
    
    %Just load the computed coefficients
    alpha1=0; 
    alpha2=-0.013208894074369; 
    alpha3=0;
    alpha4=-1.817206413939849e-04; 
    alpha5=0.001339804485786;
    alpha6=6.071198284057354e-06;
    
    beta1=2.924930866675576e-04;
    beta2=-0.011033010870322;
    beta3=-5.634557791619842e-07;
    beta4=-2.347404625178936e-05;
    beta5=0.001146515439908;
    beta6=6.906004360487817e-06;
    
    gamma1=0;
    gamma2=0.012843280276767;
    gamma3=0;
    gamma4=-0.001339100393761;
    gamma5=-1.698608716476552e-04;
else
% compute material parameters
K=@(y)(K1+K2)/2+abs(K1-K2)/2*sin(2*pi*y);
p=@(y)1./K(y); 
%p=@(y)(p1+p2)/2+abs(p1-p2)/2*sin(2*pi*y);

pm=quad(p,0,1)
ph=quad(@(y)1./p(y),0,1)^-1
Kh=quad(@(y)1./K(y),0,1)^-1

disp('computing A')
Ap=@(y,A) Kh./K(y)-ph./p(y);
[y,A]=bvp(Ap);
plot(y,A); pause(1)

disp('computing B')
Bp=@(y,B) Kh./K(y)-1;
[~,B]=bvp(Bp);
plot(y,B);  pause(1)

disp('computing C')
Cp=@(y,C) p(y)*pm^-1-1;
[~,C]=bvp(Cp);
plot(y,C); pause(1)

disp('computing D')
Dp=@(yd,D) Kh./K(yd).*get(y,C,yd)-ph./p(yd).*get(y,C,yd)-get(y,A,yd);
[~,D]=bvp(Dp);
plot(y,D); pause(1)

disp('computing E')
Ep=@(yd,E) Kh./K(yd).*get(y,C,yd)-get(y,B,yd);
[~,E]=bvp(Ep);
plot(y,E); pause(1)

disp('computing F')
Fp=@(yd,F) pm^-1*p(yd).*get(y,B,yd)-get(y,C,yd);
[~,F]=bvp(Fp);
plot(y,F); pause(1)

disp('computing H')
Hp=@(yd,H) ph^-1*p(yd).*get(y,A,yd);
[~,H]=bvp(Hp);
plot(y,H); pause(1)

disp('computing I')
auxK=average(y,K(y).^-1.*F);
auxp=average(y,p(y).^-1.*F);
Ip=@(yd,I) K(yd).^-1*Kh.*(get(y,F,yd)-Kh*auxK)...
                    -get(y,D,yd)...
                    -p(yd).^-1*ph*(get(y,F,yd)-ph*auxp);
[~,I]=bvp(Ip);
plot(y,I); pause(1)

disp('computing J')
auxK=average(y,K(y).^-1.*H);
auxp=average(y,p(y).^-1.*H);
Jp=@(yd,J) K(yd).^-1*Kh.*(get(y,H,yd)-Kh*auxK)...
                    -p(yd).^-1*ph*(get(y,H,yd)-ph*auxp);
[~,J]=bvp(Jp);
plot(y,J); pause(1)

disp('computing L')
auxK=average(y,K(y).^-1.*H);
Lp=@(yd,L) K(yd).^-1*Kh.*(get(y,H,yd)-Kh*auxK);
[~,L]=bvp(Lp);
plot(y,L); pause(1)

disp('computing M')
auxK=average(y,K(y).^-1.*F);
Mp=@(yd,M) K(yd).^-1*Kh.*(get(y,F,yd)-Kh*auxK)-get(y,E,yd);
[~,M]=bvp(Mp);
plot(y,M); pause(1)

disp('computing N')
auxp=average(y,p(y).*E);
Np=@(yd,N) p(yd)*pm^-1.*(get(y,E,yd)-pm^-1*auxp)-get(y,F,yd);
[~,N]=bvp(Np);
plot(y,N); pause(1)

disp('computing P')
auxp=average(y,p(y).*D);
Pp=@(yd,P) p(yd)*ph^-1.*(get(y,D,yd)-pm^-1*auxp)-get(y,H,yd);
[~,P]=bvp(Pp);
plot(y,P); pause(1)

disp('computing Q')
auxK=average(y,K(y).^-1.*F);
auxp=average(y,p(y).^-1.*F);
Qp=@(yd,Q) K(yd).^-1*Kh.*(get(y,N,yd)-Kh*get(y,C,yd)*auxK)...
                    -p(yd).^-1*ph.*(get(y,N,yd)-ph*get(y,C,yd)*auxp)...
                    -get(y,I,yd);
[~,Q]=bvp(Qp);
plot(y,Q); pause(1)

disp('computing R')
auxK=average(y,K(y).^-1.*H);
auxp=average(y,p(y).^-1.*H);
Rp=@(yd,R) K(yd).^-1*Kh.*(get(y,P,yd)-Kh*get(y,C,yd)*auxK)...
                    -p(yd).^-1*ph.*(get(y,P,yd)-ph*get(y,C,yd)*auxp)...
                    -get(y,J,yd);
[~,R]=bvp(Rp);
plot(y,R); pause(1)

disp('computing S')
auxK=average(y,K(y).^-1.*H);
Sp=@(yd,S) K(yd).^-1*Kh.*(get(y,P,yd)-Kh*get(y,C,yd)*auxK)...
                    -get(y,L,yd);
[~,S]=bvp(Sp);
plot(y,S); pause(1)

disp('computing T')
auxK=average(y,K(y).^-1.*F);
Tp=@(yd,S) K(yd).^-1*Kh.*(get(y,N,yd)-Kh*get(y,C,yd)*auxK)...
                    -get(y,M,yd);
[~,T]=bvp(Tp);
plot(y,T); pause(1)

disp('computing U')
auxp=average(y,p(y).*E);
Up=@(yd,U) p(yd)*pm^-1.*(get(y,M,yd)-pm^-1*get(y,B,yd)*auxp)...
                    -get(y,N,yd);
[~,U]=bvp(Up);
plot(y,U); pause(1)

disp('computing V')
auxp1=average(y,p(y).^-1.*F);
auxp2=average(y,p(y).*D);
Vp=@(yd,V) auxp1*p(yd).*get(y,A,yd)+ph^-1*p(yd).*get(y,I,yd)...
                    +pm^-1*p(yd).*(get(y,L,yd)-ph^-1*get(y,B,yd)*auxp2)...
                    -get(y,P,yd);
[~,V]=bvp(Vp);
plot(y,V); pause(1)

disp('computing W')
auxp=average(y,p(y).^-1.*H);
Wp=@(yd,W) auxp*p(yd).*get(y,A,yd)+ph^-1*p(yd).*get(y,J,yd);
[~,W]=bvp(Wp);
plot(y,W); pause(1)

disp('computing At')
auxK1=average(y,K(y).^-1.*W);
auxK2=average(y,K(y).^-1.*H);
auxp1=average(y,p(y).^-1.*W);
auxp2=average(y,p(y).^-1.*H);
Atp=@(yd,At) K(yd).^-1*Kh.*(get(y,W,yd)-Kh*auxK1-get(y,H,yd)*Kh*auxK2+Kh^2*auxK2^2)...
                -p(yd).^-1*ph.*(get(y,W,yd)-ph*auxp1-get(y,H,yd)*ph*auxp2+ph^2*auxp2^2);
[~,At]=bvp(Atp);
plot(y,At); pause(1)

disp('computing Bt')
auxp1=average(y,p(y).^-1.*W);
auxp2=average(y,p(y).^-1.*H);
Btp=@(yd,Bt) auxp1*p(yd).*get(y,A,yd)...
                    +auxp2*p(yd).*get(y,J,yd)+ph^-1*p(yd).*get(y,At,yd);
[~,Bt]=bvp(Btp);
plot(y,Bt); pause(1)

% compute coefficients
disp('*********************************')
disp('*********************************')
disp('****** COMPUTE COEFFICIENTS *****')
disp('*********************************')
disp('*********************************')
disp('')
% compute first correction coefficients
disp('COEFFICIENTS FOR FIRST CORRECTION')
[average(y,K(y).^-1.*C) average(y,p(y).^-1.*C) average(y,p(y).*A) average(y,p(y).*B)]

% compute third correction coefficients
disp('COEFFICIENTS FOR THIRD CORRECTION')
[average(y,K(y).^-1.*C) average(y,K(y).^-1.*N) average(y,K(y).^-1.*P) average(y,p(y).^-1.*C)...
    average(y,p(y).^-1.*N) average(y,p(y).^-1.*P) average(y,p(y).*A) average(y,p(y).*B)...
    average(y,p(y).*I) average(y,p(y).*J) average(y,p(y).*L) average(y,p(y).*M)]

% compute fifth correction coefficients
[average(y,p(y).*A) average(y,p(y).*J) average(y,p(y).*At)]

format long
% c-dispersion
alpha1=Kh*average(y,K(y).^-1.*F)
alpha2=Kh*average(y,K(y).^-1.*H)
alpha3=Kh*average(y,K(y).^-1.*U)-Kh^2*average(y,K(y).^-1.*F)^2
alpha4=Kh*average(y,K(y).^-1.*W)-Kh^2*average(y,K(y).^-1.*H)^2
alpha5=Kh*average(y,K(y).^-1.*V)-2*Kh^2*average(y,K(y).^-1.*F)*average(y,K(y).^-1.*H)
alpha6=Kh*average(y,K(y).^-1.*Bt)+Kh^3*average(y,K(y).^-1.*H)^3-2*Kh^2*average(y,K(y).^-1.*H)*average(y,K(y).^-1.*W)

beta1=-ph*average(y,p(y).^-1.*F)
beta2=-ph*average(y,p(y).^-1.*H)
beta3=-ph*average(y,p(y).^-1.*U)
beta4=-ph*average(y,p(y).^-1.*W)
beta5=-ph*average(y,p(y).^-1.*V)
beta6=-ph*average(y,p(y).^-1.*Bt)

gamma1=pm^-1*average(y,p(y).*E)
gamma2=ph^-1*average(y,p(y).*D)
gamma3=pm^-1*average(y,p(y).*T)-pm^-2*average(y,p(y).*E)^2
gamma4=ph^-1*average(y,p(y).*Q)+average(y,p(y).^-1.*F)*average(y,p(y).*D)+pm^-1*average(y,p(y).*S)-pm^-1*ph^-1*average(y,p(y).*E)*average(y,p(y).*D)
gamma5=average(y,p(y).*D)*average(y,p(y).^-1.*H)+ph^-1*average(y,p(y).*R)

end

function [x,y] = bvp(fun)
    x=linspace(0,1,100);
    init.x=x; init.y=0*x;
    bc=@(l,r) l(1)-1;
    sol = bvp4c(fun,bc,init);
    mean=average(sol.x,sol.y(1,:));
    x=sol.x; y=sol.y(1,:)-mean;
end

function average = average(y,F)
    average=quad(@(yd)get(y,F,yd),0,1,eps);
end

function interp = get(y,F,yi)
    interp = interp1(y,F,yi,'cubic');
end


end