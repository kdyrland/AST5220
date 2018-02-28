%% Milestone 1
% Karianne Dyrland

load cmbparams.mat

%% params

% Output options
lmax = 1200;

% Units
eV  = 1.60217646e-19;
Mpc = 3.08568025e22;
Gpc = Mpc*1e3;

% Cosmological parameters now
omega_b = 0.046;        % baryons
omega_m = 0.224;        % matter
omega_r = 8.3e-5;       % radiation
omega_nu = 0;           % neutrinos - not included in this project
omega_v = 1 - omega_m - omega_b - omega_r - omega_nu; % dark matter

T0 = 2.725;
n_s = 1;
A_s = 1;
h0 = 0.7;
H0 = h0 * 100 * 1e3 / Mpc;

% General constants
c = 2.99792458e8;               % [m/s^2]
eps_0 = 13.605698 * eV;
m_e = 9.10938188e-31;
m_H = 1.673534e-27;
sigma_T = 6.652462e-29;
G_grav = 6.67258e-11;
rho_c = 3*H0^2 / (8*pi*G_grav);
alpha = 7.29735308e-3;
hbar = 1.05457148e-34;
k_b = 1.3806503e-23;

% denisty parameters today
rho_m0 = omega_m * rho_c;
rho_b0 = omega_b * rho_c;
rho_r0 = omega_r * rho_c;
rho_v0 = omega_v * rho_c;

%% time_mod

n1 = 200;   % grid values
n2 = 300;

% redshift z
zri = 1630.4;   % start recombination
zrf = 614.2;    % end recombination
z0 = 0;         % today

% x values
xri = -log(1 + zri);    % x = ln a = -ln (1 + z)
xrf = -log(1 + zrf);
x0 = -log(1 + z0);

xt_rec = linspace(xri, xrf, n1);    % during rec
xt_today = linspace(xrf, x0, n2);   % from rec to today

xt = [xt_rec xt_today];             % concentrated array

% scale factor a
ari = 1/(1 + zri);
arf = 1/(1 + zrf);
a0 = 1;

at_rec = linspace(ari, arf, n1);    % during rec
at_today = linspace(arf, a0, n2);   % from rec to today

at = [at_rec at_today];             % concentrated array

% x values for evaluating conformal time
neta = 1000;
a1 = 1e-10;
a2 = 1;
aeta = linspace(a1,a2,neta);

xeta1 = log(a1);
xeta2 = log(a2);
xeta = linspace(xeta1, xeta2, neta);


%% Defining Hubble
% H
H = getH(xeta);
% H prime
Hp = getHp(xeta);
% dH'dx
dHp = getdHp(xeta);


%% get eta

% derivatives of eta
dnda = etaderivs(xeta);

% integrate ode using ode45
[x,eta] = ode45(@(xeta,eta) etaderivs(xeta), xeta, xeta1);

% splined eta
etanew = spline(xeta,eta,xeta);

figure
semilogy(xeta,eta/Gpc,'Displayname','ODE 45')
hold on
semilogy(xeta,etanew/Gpc,'Displayname','Spline')
title('Conformal time','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('\eta [Mpc]','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;

%% plot hubble

figure
semilogy(xeta,H.*Mpc./1e3)
title('Hubble parameter','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('H [km/s/Mpc]','Fontsize',14)

%% plot hubble prime

zuni = linspace(zri,z0,1000);       % uniform redshift array

figure
semilogy(zuni,Hp.*Gpc.*1e3)
set(gca,'Xdir','reverse')
title('Scaled Hubble parameter','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('$\mathcal{H}$ [km/s/Mpc]','Interpreter','latex','Fontsize',14)

%% find omegas

rho_c = (3.*H.^2)./(8.*pi.*G_grav);     % critical density for all time

Omega_m = rho_m0 * exp(xeta).^(-3) ./ rho_c;
Omega_b = rho_b0 * exp(xeta).^(-3) ./ rho_c;
Omega_r = rho_r0 * exp(xeta).^(-4) ./ rho_c;
Omega_v = rho_v0 ./ rho_c;
Omega_tot = Omega_m + Omega_b + Omega_r + Omega_v;

figure
hold on
plot(xeta,Omega_b,'Displayname','\Omega_b')
plot(xeta,Omega_m,'Displayname','\Omega_m')
plot(xeta,Omega_r,'Displayname','\Omega_r')
plot(xeta,Omega_v,'Displayname','\Omega_\Lambda')
% plot(xeta,Omega_tot,'Displayname','\Omega_{total}')
title('Density parameters','Fontsize',14)
xlabel('x','Fontsize',14)
ylabel('\Omega various','Fontsize',14)

lgd = legend('show','Location','Northwest');
lgd.FontSize = 14;








