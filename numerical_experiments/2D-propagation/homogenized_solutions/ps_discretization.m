function Du = ps_discretization(u,kx,ky)

global alpha1 alpha2 alpha3 alpha4 alpha5 beta1 beta2 beta3 beta4 beta5 gamma1 gamma2 gamma3 gamma4 gamma5 Kh ph pm
delta2=1;
delta4=1;

[kx,ky]=meshgrid(kx,ky);

uhat=fft2(squeeze(u(1,:,:)));
vhat=fft2(squeeze(u(2,:,:)));
phat=fft2(squeeze(u(3,:,:)));

dudx        = real(ifft2(1i *kx             .*uhat));
du3dx3      = real(ifft2(-1i*kx.^3          .*uhat));
du3dxdy2    = real(ifft2(-1i*kx.*ky.^2      .*uhat));
du5dxdy4    = real(ifft2(1i *kx.*ky.^4      .*uhat));
du5dx5      = real(ifft2(1i *kx.^5          .*uhat));
du5dx3dy2   = real(ifft2(1i *kx.^3.*ky.^2   .*uhat));

dvdy        = real(ifft2(1i *ky              .*vhat));
dv3dy3      = real(ifft2(-1i*ky.^3          .*vhat));
dv3dx2dy    = real(ifft2(-1i*ky.*kx.^2      .*vhat));
dv5dy5      = real(ifft2(1i *ky.^5           .*vhat));
dv5dx4dy    = real(ifft2(1i *kx.^4.*ky       .*vhat));
dv5dx2dy3   = real(ifft2(1i *kx.^2.*ky.^3    .*vhat));

dpdx        = real(ifft2(1i *kx             .*phat));
dpdy        = real(ifft2(1i *ky             .*phat));
dp3dx3      = real(ifft2(-1i*kx.^3          .*phat));
dp3dy3      = real(ifft2(-1i*ky.^3          .*phat));
dp3dxdy2    = real(ifft2(-1i*kx.*ky.^2      .*phat));
dp3dx2dy    = real(ifft2(-1i*ky.*kx.^2      .*phat));
dp5dxdy4    = real(ifft2(1i *kx.*ky.^4      .*phat));
dp5dx5      = real(ifft2(1i *kx.^5          .*phat));
dp5dx3dy2   = real(ifft2(1i *kx.^3.*ky.^2   .*phat));
dp5dy5      = real(ifft2(1i *ky.^5          .*phat));
dp5dx2dy3   = real(ifft2(1i *kx.^2.*ky.^3   .*phat));
dp5dx4dy    = real(ifft2(1i *kx.^4.*ky      .*phat));

Du(1,:,:)=ph^-1*    ( -dpdx  + delta2*(beta1*dp3dxdy2+beta2*dp3dx3)... 
                             + delta4*(beta3*dp5dxdy4+beta4*dp5dx5+beta5*dp5dx3dy2) );

Du(2,:,:)=pm^-1*    ( -dpdy  + delta2*(gamma1*dp3dy3+gamma2*dp3dx2dy)...
                             + delta4*(gamma3*dp5dy5+gamma4*dp5dx2dy3+gamma5*dp5dx4dy) );
                         
Du(3,:,:)=Kh*       (-(dudx+dvdy) ...
                            + delta2*(alpha1*(du3dxdy2+dv3dy3)+alpha2*(du3dx3+dv3dx2dy))...
                            + delta4*(alpha3*(du5dxdy4+dv5dy5)+alpha4*(du5dx5+dv5dx4dy)...
                                                +alpha5*(du5dx3dy2+dv5dx2dy3)) );
                        
                        
                        
                        