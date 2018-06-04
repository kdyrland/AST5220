%% plot data

load ms4data.mat        % data imported and saved to .mat workspace

Mpc = 3.08568025e22;
h0 = 0.7;
H0 = h0 * 100.d0 * 1.d3 / Mpc;
c = 2.99792458e8;
hc = H0/c;
ch = c/H0;

%% test bessel
xh = xhires(2501:end);
int = integrand(2501:end)*1e3;

figure
plot(xh,int)
ylim([-3,2])
xlim([-8 0])
xlabel('$x$','Interpreter','latex','Fontsize',16)
ylabel('$\hat{S}(k,x) j_l[k(\eta_0 - \eta(x))]/10^{-3}$','Interpreter','latex','Fontsize',16)
lgd = legend('show','k = 340 H_0/c,  l = 100','Location','Southeast');
lgd.FontSize = 12;

figure
plot(xh,int)
ylim([-0.7,0.7])
xlim([-1.5 0])
xlabel('$x$','Interpreter','latex','Fontsize',16)
ylabel('$\hat{S}(k,x) j_l[k(\eta_0 - \eta(x))]/10^{-3}$','Interpreter','latex','Fontsize',16)


%% integrands

lw = [6 100 200 500 800 1200];
kch = khires*ch;
int1 = cl_int1.*lw(1)*(lw(1)+1)*hc;
int2 = cl_int2.*lw(2)*(lw(2)+1)*hc;
int3 = cl_int3.*lw(3)*(lw(3)+1)*hc;
int4 = cl_int4.*lw(4)*(lw(4)+1)*hc;
int5 = cl_int5.*lw(5)*(lw(5)+1)*hc;
int6 = cl_int6.*lw(6)*(lw(6)+1)*hc;

%plot
figure
hold on

plot(kch,int1,'Displayname','l = 6')
plot(kch,int2,'Displayname','l = 100')
plot(kch,int3,'Displayname','l = 200')
plot(kch,int4,'Displayname','l = 500')
plot(kch,int5,'Displayname','l = 800')
plot(kch,int6,'Displayname','l = 1200')

xlim([0 500])
lgd = legend('show');
lgd.FontSize = 14;
xlabel('$ck/H_0$','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)\Theta^2(k) H_0/ck$','Interpreter','latex','Fontsize',16)


%% cmb

lplanck = horzcat(lplanck1.',lplanck2.');
clplanck = horzcat(clplanck1.',clplanck2.');
planckerror = horzcat(planckerror1.',planckerror2.');
cl = cls.*5775./max(cls);   % normalized to planck

figure
hold on
errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
plot(lhires,cl,'Linewidth',1.2,'Displayname','Simulated','Color',[255/255 153/255 51/255])
xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 14;

%% parameter estimation - baryons

clb0 = cl_b0.*5775./max(cl_b0);
clb02 = cl_b02.*5775./max(cl_b02);
clb04 = cl_b04.*5775./max(cl_b04);
clb1 = cl_b1.*5775./max(cl_b1);

figure
hold on

errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
%plot(lhires,clb0,'Linewidth',1,'Displayname','\Omega_b = 0.0')
plot(lhires,clb02,'Linewidth',1,'Displayname','\Omega_b = 0.020')
plot(lhires,cl,'Linewidth',1,'Displayname','\Omega_b = 0.042')
plot(lhires,clb1,'Linewidth',1,'Displayname','\Omega_b = 0.100')
xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 12;

%% parameter estimation - dark matter

clb0 = cl_b0.*5775./max(cl_b0);
clm1 = cl_m1.*5775./max(cl_m1);
clm2 = cl_m2.*5775./max(cl_m2);
clm3 = cl_m3.*5775./max(cl_m3);

figure
hold on

errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
plot(lhires,clm1,'Linewidth',1,'Displayname','\Omega_m = 0.100')
plot(lhires,clm2,'Linewidth',1,'Displayname','\Omega_m = 0.200')
plot(lhires,cl,'Linewidth',1,'Displayname','\Omega_m = 0.224')
plot(lhires,clm3,'Linewidth',1,'Displayname','\Omega_m = 0.300')

xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 12;



%% parameter estimation - hubble

clh6 = cl_h6.*5775./max(cl_h6);
clh66 = cl_h66.*5775./max(cl_h66);
clh8 = cl_h08.*5775./max(cl_h08);

figure
hold on

errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
plot(lhires,clh6,'Linewidth',1,'Displayname','h = 0.60')
plot(lhires,clh66,'Linewidth',1,'Displayname','h = 0.66')
plot(lhires,cl,'Linewidth',1,'Displayname','h = 0.70')
plot(lhires,clh8,'Linewidth',1,'Displayname','h = 0.80')

xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 12;


%% parameter estimation - n_s

cln06 = cl_n06.*5775./max(cl_n06);
cln07 = cl_n07.*5775./max(cl_n07);
cln08 = cl_n08.*5775./max(cl_n08);
cln12 = cl_n12.*5775./max(cl_n12);

figure
hold on

errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
plot(lhires,cln07,'Linewidth',1,'Displayname','n_s = 0.7')
plot(lhires,cln08,'Linewidth',1,'Displayname','n_s = 0.8')
plot(lhires,cl,'Linewidth',1,'Displayname','n_s = 1.0')
plot(lhires,cln12,'Linewidth',1,'Displayname','n_s = 1.2')

xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 12;


%% best fit

clfit = cl_fit.*5775./max(cl_fit);

figure
hold on
errorbar(lplanck,clplanck,planckerror,'Color',[102/255 153/255 153/255],'Displayname','Planck data')
plot(lhires,cl,'Linewidth',1.2,'Displayname','Simulated','Color',[255/255 153/255 51/255])
plot(lhires,clfit,'Linewidth',1,'Displayname','Best fit')
xlim([0 1200])
xlabel('\it{l}','Interpreter','latex','Fontsize',16)
ylabel('$l(l+1)C_l/2\pi$','Interpreter','latex','Fontsize',16)

lgd = legend('show');
lgd.FontSize = 14;
















