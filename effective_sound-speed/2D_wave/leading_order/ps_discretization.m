function Du = ps_discretization(u,kx,ky)

global Kh ph pm

[kx,ky]=meshgrid(kx,ky);

uhat=fft2(squeeze(u(1,:,:)));
vhat=fft2(squeeze(u(2,:,:)));
phat=fft2(squeeze(u(3,:,:)));

dudx        = real(ifft2(1i *kx             .*uhat));
dvdy        = real(ifft2(1i *ky              .*vhat));

dpdx        = real(ifft2(1i *kx             .*phat));
dpdy        = real(ifft2(1i *ky             .*phat));

Du(1,:,:)=ph^-1*    ( -dpdx );

Du(2,:,:)=pm^-1*    ( -dpdy );
                         
Du(3,:,:)=Kh*       (-(dudx+dvdy) );
                        
                        
                        
                        