%% plot data

load ms2data.mat        % data imported and saved to .mat workspace

%% electron density
arec = exp(xrec);
zrec = 1./arec - 1;

semilogy(zrec,Xe)
xlim([0 1800])
set(gca,'XDir','reverse')
xlabel('z','FontSize',14)
ylabel('X_e(z)','FontSize',14)
title('Fractional electron density','FontSize',14)


%% optical depth

figure
semilogy(tau,'Displayname','\tau')
hold on
semilogy(abs(dtau),'Displayname','\tau''')
semilogy(tau2, 'Displayname','\tau''''')
xlabel('x','FontSize',14);
ylabel('Optical Depth','FontSize',14)
lgd = legend('show');
lgd.FontSize = 14;


%% visibility function

figure
plot(g/max(g),'Displayname','g')
hold on
plot(g2/max(g2),'Displayname','g''')
plot(dg/max(dg),'Displayname','g''''')
xlim([660 740])
ylim([-1.5 1.2])
xlabel('x','FontSize',14)
ylabel('Visibility function','FontSize',14)
lgd = legend('show');
lgd.FontSize = 14;










