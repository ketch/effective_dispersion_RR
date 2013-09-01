function parameters(KA,pA,KB,pB)

global alpha1 alpha2 alpha3 alpha4 alpha5 beta1 beta2 beta3 beta4 beta5 gamma1 gamma2 gamma3 gamma4 gamma5 Kh ph pm

%ZA=sqrt(KA*pA); ZB=sqrt(KB*pB);
%cA=sqrt(KA/pA); cB=sqrt(KB/pB);

%Km=(KA+KB)/2;
Kh=2*KA*KB/(KA+KB)
pm=(pA+pB)/2
ph=2*pA*pB/(pA+pB)

%alpha1=-(KA-KB)/(192*Km^2)*(ZA^2-ZB^2)/pm;
%alpha2=-(KA-KB)/(192*Km^2)*(cA^2-cB^2)*pm;

%beta1=(pA-pB)/(192*Km)*(ZA^2-ZB^2)/pm^2;
%beta2=(pA-pB)/(192*Km)*(cA^2-cB^2);

%gamma1=-(pA-pB)/(192*Km)*(ZA^2-ZB^2)/pm^2;
%gamma2=-(pA-pB)/(192*Km)*(cA^2-cB^2);

alpha1=(-1/24).*(KA+(-1).*KB).*(KA+KB).^(-2).*(pA+pB).^(-1).*(KA.*pA+(-1).*KB.*pB);
alpha2=(-1/96).*(KA+(-1).*KB).*(KA+KB).^(-2).*pA.^(-1).*pB.^(-1).*(pA+pB).*((-1).*KB.*pA+KA.*pB);
alpha3=(-1/5760).*(KA+(-1).*KB).*(KA+KB).^(-4).*(pA+pB).^(-3).*(KA.*pA+(-1).*KB.*pB).*(KA.^2.*(3.*pA.^2+(-38).*pA.*pB+(-21).*pB.^2)+...
  KB.^2.*((-21).*pA.^2+(-38).*pA.*pB+3.*pB.^2)+(-2).*KA.*KB.*(29.*pA.^2+78.*pA.*pB+29.*pB.^2));
alpha4=(-1/11520).*(KA+(-1).*KB).*(KA+KB).^(-4).*pA.^(-2).*pB.^(-2).*(pA+...
    pB).*(KB.*pA+(-1).*KA.*pB).^2.*(3.*KA.*pA+(-2).*KB.*pA+2.*KA.*pB+(-3).*KB.*pB);
alpha5=(1/552960).*((-1).*KA+KB).*(KA+KB).^(-5).*pA.^(-1).*pB.^(-1).*(pA+...
    pB).^(-1).*(KB.*pA+(-1).*KA.*pB).*(24.*KB.^3.*(21.*pA.^2+22.*pA.*pB+(-27).*pB.^2)+KA.^3.*((-147).*pA.^2+1010.*pA.*pB+505.*pB.^2)+ ...
    KA.*KB.^2.*(2281.*pA.^2+4938.*pA.*pB+677.*pB.^2)+2.*KA.^2.*KB.*(805.*pA.^2+2690.*pA.*pB+897.*pB.^2));

%alpha6=(-1/30965760).*(KA+(-1).*KB).*(KA+KB).^(-6).*pA.^(-3).*pB.^(-3).*(...
%    (-1).*KB.*pA+KA.*pB).^3.*(7.*KB.*pA+(204.*KA.^2+(-383).*KA.*KB+43.*KB.^2).*pA.^3+7.*(KB+(68.*KA.^2+(-167).*KA.*KB+35.*KB.^2).* ...
%    pA.^2).*pB+7.*(50.*KA.^2+(-171).*KA.*KB+49.*KB.^2).*pA.*pB.^2+7.*(10.*KA.^2+(-61).*KA.*KB+19.*KB.^2).*pB.^3);

beta1=(1/24).*(KA+KB).^(-1).*(pA+(-1).*pB).*(pA+pB).^(-2).*(KA.*pA+(-1).*KB.*pB);
beta2=(-1/96).*(KA+KB).^(-1).*pA.^(-1).*(pA+(-1).*pB).*pB.^(-1).*(KB.*pA+(-1).*KA.*pB);
beta3=(-1/5760).*(KA+KB).^(-3).*(pA+(-1).*pB).*(pA+pB).^(-4).*(KA.*pA+(-1).*KB.*pB).*(8.*KA.*KB.*(6.*pA.^2+17.*pA.*pB+6.*pB.^2)+KB.^2.*( ...
    21.*pA.^2+48.*pA.*pB+7.*pB.^2)+KA.^2.*(7.*pA.^2+48.*pA.*pB+21.*pB.^2));
beta4=(-1/46080).*(KA+KB).^(-3).*pA.^(-2).*(pA+(-1).*pB).*pB.^(-2).*(...
    KB.*pA+(-1).*KA.*pB).^2.*((-7).*KA.*pA+3.*KB.*pA+(-3).*KA.*pB+7.*KB.*pB);
beta5=(1/23040).*(KA+KB).^(-3).*pA.^(-1).*(pA+(-1).*pB).*pB.^(-1).*(pA+pB).^(-2).*(KB.*pA+(-1).*KA.*pB).*((-7).*KA.^2.*(pA.^2+(-6).*pA.* ...
    pB+(-3).*pB.^2)+7.*KB.^2.*(3.*pA.^2+6.*pA.*pB+(-1).*pB.^2)+2.*KA.*KB.*(27.*pA.^2+82.*pA.*pB+27.*pB.^2));

%beta6=(-1/30965760).*(KA+KB).^(-5).*pA.^(-3).*pB.^(-3).*((-1).*pA+pB).*( ...
%    pA+pB).^(-1).*(KB.*pA+(-1).*KA.*pB).^3.*((-7).*KB.*pA+((-71).*KA.^2+173.*KA.*KB+34.*KB.^2).*pA.^3+(-7).*(KB+(19.*KA.^2+(-77).* ...
%    KA.*KB+(-6).*KB.^2).*pA.^2).*pB+(-63).*KA.*(KA+(-9).*KB).*pA.*pB.^2+7.*KA.*(KA+31.*KB).*pB.^3);

gamma1=(-1/24).*(KA+KB).^(-1).*(pA+(-1).*pB).*(pA+pB).^(-2).*(KA.*pA+(-1).*KB.*pB);
gamma2=(-1/96).*(KA+KB).^(-1).*pA.^(-1).*(pA+(-1).*pB).*pB.^(-1).*((-1).*KB.*pA+KA.*pB);
gamma3=(-1/5760).*(KA+KB).^(-3).*(pA+(-1).*pB).*(pA+pB).^(-4).*(KA.*pA+(-1).*KB.*pB).*(KA.^2.*(3.*pA.^2+(-58).*pA.*pB+(-21).*pB.^2)+ ...
    KB.^2.*((-21).*pA.^2+(-58).*pA.*pB+3.*pB.^2)+(-2).*KA.*KB.*(19.*pA.^2+78.*pA.*pB+19.*pB.^2));
gamma4=(-1/23040).*(KA+KB).^(-3).*pA.^(-1).*(pA+(-1).*pB).*pB.^(-1).*(pA+pB).^(-2).*(KB.*pA+(-1).*KA.*pB).*(KB.^2.*(21.*pA.^2+62.*pA.*pB+( ...
    -27).*pB.^2)+34.*KA.*KB.*(pA.^2+6.*pA.*pB+pB.^2)+KA.^2.*((-27).*pA.^2+62.*pA.*pB+21.*pB.^2));
gamma5=(-1/23040).*(KA+KB).^(-3).*pA.^(-2).*(pA+(-1).*pB).*pB.^(-2).*(KB.*pA+(-1).*KA.*pB).^2.*((6.*KA+KB).*pA+(-1).*(KA+6.*KB).*pB);


