function parameters(KA,pA,KB,pB)

global Kh ph pm

%ZA=sqrt(KA*pA); ZB=sqrt(KB*pB);
%cA=sqrt(KA/pA); cB=sqrt(KB/pB);

%Km=(KA+KB)/2;
Kh=2*KA*KB/(KA+KB)
pm=(pA+pB)/2
ph=2*pA*pB/(pA+pB)

