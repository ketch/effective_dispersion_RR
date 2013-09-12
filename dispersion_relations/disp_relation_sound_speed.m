function disp_relation
clc

K1=1; rho1=1; K2=1; rho2=1;
get_plots(K1,rho1,K2,rho2,'homogeneous')

K1=5/8; rho1=8/5; K2=5/2; rho2=2/5;
get_plots(K1,rho1,K2,rho2,'c-medium')

K1=1; rho1=1; K2=4; rho2=4;
get_plots(K1,rho1,K2,rho2,'z-medium')

K1=1; rho1=1; K2=16; rho2=1;
get_plots(K1,rho1,K2,rho2,'isotropic')

rc=2.15; rz=8;
K1=1; rho1=1;
K2=rc/rz; rho2=1/(rz*rc);
get_plots(K1,rho1,K2,rho2,'cz-medium')


end

function get_plots(K1,rho1,K2,rho2,name)
global d pm ph Kh a1 a2 a3 a4 a5 a6 b1 b2 b3 b4 b5 b6 g1 g2 g3 g4 g5

parameters(K1,rho1,K2,rho2)
%disp relation
w2=@(kkx,kky)(ph.^(-1).*pm.^(-1).*(Kh.*kky.^2.*ph+Kh.*kkx.^2.*pm)...
    +d.^2.*ph.^(-1).*pm.^(-1).*(a2.*Kh.*kkx.^2.*kky.^2.*ph+g2.*Kh.*kkx.^2.*kky.^2.* ...
        ph+a1.*Kh.*kky.^4.*ph+g1.*Kh.*kky.^4.*ph+a2.*Kh.*kkx.^4.*pm+b2.*Kh.*kkx.^4.*pm+a1.*Kh.*kkx.^2.*kky.^2.*pm+b1.*Kh.*kkx.^2.*kky.^2.*pm)...
    +1*d.^4.*ph.^(-1).*pm.^(-1).*(a4.*Kh.*kkx.^4.*kky.^2.*ph+a2.*g2.*Kh.*kkx.^4.*kky.^2.*ph+(-1).*g5.*Kh.*kkx.^4.*kky.^2.*ph+(-1).*a5.* ...
        Kh.*kkx.^2.*kky.^4.*ph+a2.*g1.*Kh.*kkx.^2.*kky.^4.*ph+a1.*g2.*Kh.*kkx.^2.*kky.^4.*ph+(-1).*g4.*Kh.*kkx.^2.*kky.^4.*ph+a3.*Kh.* ...
        kky.^6.*ph+a1.*g1.*Kh.*kky.^6.*ph+(-1).*g3.*Kh.*kky.^6.*ph+a4.*Kh.*kkx.^6.*pm+a2.*b2.*Kh.*kkx.^6.*pm+(-1).*b4.*Kh.*kkx.^6.*pm+( ...
        -1).*a5.*Kh.*kkx.^4.*kky.^2.*pm+a2.*b1.*Kh.*kkx.^4.*kky.^2.*pm+a1.*b2.*Kh.*kkx.^4.*kky.^2.*pm+(-1).*b5.*Kh.*kkx.^4.*kky.^2.*pm+ ...
        a3.*Kh.*kkx.^2.*kky.^4.*pm+a1.*b1.*Kh.*kkx.^2.*kky.^4.*pm+(-1).*b3.*Kh.*kkx.^2.*kky.^4.*pm));
c=@(kkx,kky)sqrt(w2(kkx,kky))./sqrt(kkx.^2+kky.^2);

[kkx,kky]=meshgrid(linspace(0,3,500),linspace(0,3,500));
c=c(kkx,kky);

cm = cbrewer('div','RdBu',20);
colormap(cm)
colormap(flipud(colormap))

%% pcolor
f=figure(1); clf
surf(kkx,kky,c,'EdgeColor','none')
set(gca,'fontsize',25)
xlabel('k k_x'); ylabel('k k_y')
view(2); colorbar('fontsize',25)
axis equal
axis([0 3 0 3])
saveas(f,['disp-rel_speed_' name '_pcolor'],'png')

%% countour
f=figure(2); clf
contour(kkx,kky,c)
set(gca,'fontsize',25)
xlabel('k k_x'); ylabel('k k_y')
view(2); colorbar('fontsize',25)
axis equal
axis([0 3 0 3])
saveas(f,['disp-rel_speed_' name '_contour'],'png')

%% plot slices
f=figure(3); clf; hold on
plot(kkx(1,:),c(1,:),'b','linewidth',3)
plot(kky(:,1),c(:,1),'--r','linewidth',3)
set(gca,'fontsize',20)
xlabel('k k_x, k k_y'); ylabel('c=\omega/k')
legend('Transverse propagation', 'Normal propagation','location','best')
saveas(f,['disp-rel_speed_' name '_slices'],'epsc')

end
