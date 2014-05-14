function coeff
clc; clear all; clf
addpath '/Users/mquezada/Desktop/dispersion-relation/chebfun'

% compute material parameters
K1=5/8; p1=8/5;
K2=5/2; p2=2/5;
K=@(y)(K1+K2)/2+abs(K1-K2)/2*sin(2*pi*y);
p=@(y)1./K(y); 

pm=quad(p,0,1);
ph=quad(@(y)1./p(y),0,1)^-1
Kh=quad(@(y)1./K(y),0,1)^-1

disp('computing A')
Ap=@(y,A) Kh./K(y)-ph./p(y);
[y,A]=bvp(Ap);
plot(y,A); pause(1)

disp('computing B')
Bp=@(y,B) Kh./K(y)-1;
[dummy,B]=bvp(Bp);
plot(y,B);  pause(1)

disp('computing C')
Cp=@(y,C) p(y)*pm^-1-1;
[dummy,C]=bvp(Cp);
plot(y,C); pause(1)


disp('computing D')
Dp=@(yd,D) Kh./K(yd).*get(y,C,yd)-ph./p(yd).*get(y,C,yd)-get(y,A,yd);
[dummy,D]=bvp(Dp);
plot(y,D); pause(1)

disp('computing E')
Ep=@(yd,E) Kh./K(yd).*get(y,C,yd)-get(y,B,yd);
[dummy,E]=bvp(Ep);
plot(y,E); pause(1)

disp('computing F')
Fp=@(yd,F) pm^-1*p(yd).*get(y,B,yd)-get(y,C,yd);
[dummy,F]=bvp(Fp);
plot(y,F); pause(1)

disp('computing H')
Hp=@(yd,H) ph^-1*p(yd).*get(y,A,yd);
[dummy,H]=bvp(Hp);
plot(y,H); pause(1)

disp('computing I')
auxK=average(y,K(y).^-1.*F);
auxp=average(y,p(y).^-1.*F);
Ip=@(yd,I) K(yd).^-1*Kh.*(get(y,F,yd)-Kh*auxK)...
                    -get(y,D,yd)...
                    -p(yd).^-1*ph*(get(y,F,yd)-ph*auxp);
[dummy,I]=bvp(Ip);
plot(y,I); pause(1)

disp('computing J')
auxK=average(y,K(y).^-1.*H);
auxp=average(y,p(y).^-1.*H);
Jp=@(yd,J) K(yd).^-1*Kh.*(get(y,H,yd)-Kh*auxK)...
                    -p(yd).^-1*ph*(get(y,H,yd)-ph*auxp);
[dummy,J]=bvp(Jp);
plot(y,J); pause(1)

disp('computing L')
auxK=average(y,K(y).^-1.*H);
Lp=@(yd,L) K(yd).^-1*Kh.*(get(y,H,yd)-Kh*auxK);
[dummy,L]=bvp(Lp);
plot(y,L); pause(1)

disp('computing M')
auxK=average(y,K(y).^-1.*F);
Mp=@(yd,M) K(yd).^-1*Kh.*(get(y,F,yd)-Kh*auxK)-get(y,E,yd);
[dummy,M]=bvp(Mp);
plot(y,M); pause(1)

disp('computing N')
auxp=average(y,p(y).*E);
Np=@(yd,N) p(yd)*pm^-1.*(get(y,E,yd)-pm^-1*auxp)-get(y,F,yd);
[dummy,N]=bvp(Np);
plot(y,N); pause(1)

disp('computing P')
auxp=average(y,p(y).*D);
Pp=@(yd,P) p(yd)*ph^-1.*(get(y,D,yd)-pm^-1*auxp)-get(y,H,yd);
[dummy,P]=bvp(Pp);
plot(y,P); pause(1)

disp('computing Q')
auxK=average(y,K(y).^-1.*F);
auxp=average(y,p(y).^-1.*F);
Qp=@(yd,Q) K(yd).^-1*Kh.*(get(y,N,yd)-Kh*get(y,C,yd)*auxK)...
                    -p(yd).^-1*ph.*(get(y,N,yd)-ph*get(y,C,yd)*auxp)...
                    -get(y,I,yd);
[dummy,Q]=bvp(Qp);
plot(y,Q); pause(1)

disp('computing R')
auxK=average(y,K(y).^-1.*H);
auxp=average(y,p(y).^-1.*H);
Rp=@(yd,R) K(yd).^-1*Kh.*(get(y,P,yd)-Kh*get(y,C,yd)*auxK)...
                    -p(yd).^-1*ph.*(get(y,P,yd)-ph*get(y,C,yd)*auxp)...
                    -get(y,J,yd);
[dummy,R]=bvp(Rp);
plot(y,R); pause(1)

disp('computing S')
auxK=average(y,K(y).^-1.*H);
Sp=@(yd,S) K(yd).^-1*Kh.*(get(y,P,yd)-Kh*get(y,C,yd)*auxK)...
                    -get(y,L,yd);
[dummy,S]=bvp(Sp);
plot(y,S); pause(1)

disp('computing T')
auxK=average(y,K(y).^-1.*F);
Tp=@(yd,S) K(yd).^-1*Kh.*(get(y,N,yd)-Kh*get(y,C,yd)*auxK)...
                    -get(y,M,yd);
[dummy,T]=bvp(Tp);
plot(y,T); pause(1)

disp('computing U')
auxp=average(y,p(y).*E);
Up=@(yd,U) p(yd)*pm^-1.*(get(y,M,yd)-pm^-1*get(y,B,yd)*auxp)...
                    -get(y,N,yd);
[dummy,U]=bvp(Up);
plot(y,U); pause(1)

disp('computing V')
auxp1=average(y,p(y).^-1.*F);
auxp2=average(y,p(y).*D);
Vp=@(yd,V) auxp1*p(yd).*get(y,A,yd)+ph^-1*p(yd).*get(y,I,yd)...
                    +pm^-1*p(yd).*(get(y,L,yd)-ph^-1*get(y,B,yd)*auxp2)...
                    -get(y,P,yd);
[dummy,V]=bvp(Vp);
plot(y,V); pause(1)

disp('computing W')
auxp=average(y,p(y).^-1.*H);
Wp=@(yd,W) auxp*p(yd).*get(y,A,yd)+ph^-1*p(yd).*get(y,J,yd);
[dummy,W]=bvp(Wp);
plot(y,W); pause(1)

disp('computing At')
auxK1=average(y,K(y).^-1.*W);
auxK2=average(y,K(y).^-1.*H);
auxp1=average(y,p(y).^-1.*W);
auxp2=average(y,p(y).^-1.*H);
Atp=@(yd,At) K(yd).^-1*Kh.*(get(y,W,yd)-Kh*auxK1-get(y,H,yd)*Kh*auxK2+Kh^2*auxK2^2)...
                -p(yd).^-1*ph.*(get(y,W,yd)-ph*auxp1-get(y,H,yd)*ph*auxp2+ph^2*auxp2^2);
[dummy,At]=bvp(Atp);
plot(y,At); pause(1)

disp('computing Bt')
auxp1=average(y,p(y).^-1.*W);
auxp2=average(y,p(y).^-1.*H);
Btp=@(yd,Bt) auxp1*p(yd).*get(y,A,yd)...
                    +auxp2*p(yd).*get(y,J,yd)+ph^-1*p(yd).*get(y,At,yd);
[dummy,Bt]=bvp(Btp);
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
disp('COEFFICIENTS FOR FIFTH CORRECTION')
[average(y,p(y).*A) average(y,p(y).*J) average(y,p(y).*At)]

format long
% c-dispersion
alpha2=Kh*average(y,K(y).^-1.*H)
alpha4=-Kh^2*average(y,K(y).^-1.*H)^2+Kh*average(y,K(y).^-1.*W)
alpha6=-2*Kh^2*average(y,K(y).^-1.*H)*average(y,K(y).^-1.*W)+Kh^3*average(y,K(y).^-1.*H)^3 ...
            +Kh*average(y,K(y).^-1.*Bt)

beta2=-ph*average(y,p(y).^-1.*H)
beta5=-ph*average(y,p(y).^-1.*W)
beta6=-ph*average(y,p(y).^-1.*Bt)

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
