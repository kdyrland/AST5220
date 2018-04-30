%% plot data

load ms2data.mat        % data imported and saved to .mat workspace

%% electron density
arec = exp(xrec);
zrec = 1./arec - 1;

semilogy(zrec,Xe)
xlim([0 1800])
ylim([1e-4 2])
set(gca,'XDir','reverse')
xlabel('z','FontSize',14)
ylabel('X_e(z)','FontSize',14)
title('Fractional electron density','FontSize',14)


%% optical depth

figure
semilogy(xrec,tau,'Displayname','\tau')
hold on
semilogy(xrec,abs(dtau),'Displayname','\tau''')
semilogy(xrec,tau2, 'Displayname','\tau''''')
xlabel('x','FontSize',14);
ylabel('Optical Depth','FontSize',14)
lgd = legend('show');
lgd.FontSize = 14;


%% visibility function
close all
figure
hold on
plot(xrec,g,'Displayname','g')
plot(xrec,dg/10,'Displayname','g''')
plot(xrec,g2/300,'Displayname','g''''')
xlim([-7.5 -6])
% ylim([-1.5 1.2])
xlabel('x','FontSize',14)
ylabel('Visibility function','FontSize',14)
lgd = legend('show');
lgd.FontSize = 14;



