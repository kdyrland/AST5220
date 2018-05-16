%% plot data

%load ms3data.mat        % data imported and saved to .mat workspace


%% Theta_0
figure
hold on

theta01(1) = 0.5;
theta02(1) = 0.5;
theta03(1) = 0.5;
theta04(1) = 0.5;
theta05(1) = 0.5;
theta06(1) = 0.5;

plot(xf,theta01,'Displayname','k_0 = 1')
plot(xf,theta02,'Displayname','k_1 = 12')
plot(xf,theta03,'Displayname','k_3 = 30')
plot(xf,theta04,'Displayname','k_4 = 40')
plot(xf,theta05,'Displayname','k_5 = 85')
plot(xf,theta06,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\Theta_0','Fontsize',16)
lgd = legend('show','Location','Southwest');
lgd.FontSize = 12;

%% Theta_1
figure
hold on

plot(xf,theta11,'Displayname','k_0 = 1')
plot(xf,theta12,'Displayname','k_1 = 12')
plot(xf,theta13,'Displayname','k_3 = 30')
plot(xf,theta14,'Displayname','k_4 = 40')
plot(xf,theta15,'Displayname','k_5 = 85')
plot(xf,theta16,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\Theta_0','Fontsize',16)
lgd = legend('show','Location','Southwest');
lgd.FontSize = 12;

%% Phi
figure
hold on
phi1(1) = 1;
phi2(1) = 1;
phi3(1) = 1;
phi4(1) = 1;
phi5(1) = 1;
phi6(1) = 1;

plot(xf,phi1,'Displayname','k_0 = 1')
plot(xf,phi2,'Displayname','k_1 = 12')
plot(xf,phi3,'Displayname','k_3 = 30')
plot(xf,phi4,'Displayname','k_4 = 40')
plot(xf,phi5,'Displayname','k_5 = 85')
plot(xf,phi6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\Phi','Fontsize',16)
lgd = legend('show','Location','Southwest');
lgd.FontSize = 12;

%% Psi
figure
hold on

psi1(1) = -1;
psi2(1) = -1;
psi3(1) = -1;
psi4(1) = -1;
psi5(1) = -1;
psi6(1) = -1;

plot(xf,psi1,'Displayname','k_0 = 1')
plot(xf,psi2,'Displayname','k_1 = 12')
plot(xf,psi3,'Displayname','k_3 = 30')
plot(xf,psi4,'Displayname','k_4 = 40')
plot(xf,psi5,'Displayname','k_5 = 85')
plot(xf,psi6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\Psi','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;

%% v
figure
hold on
plot(xf,v1,'Displayname','k_0 = 1')
plot(xf,v2,'Displayname','k_1 = 12')
plot(xf,v3,'Displayname','k_3 = 30')
plot(xf,v4,'Displayname','k_4 = 40')
plot(xf,v5,'Displayname','k_5 = 85')
plot(xf,v6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('v','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;

%% vb
figure
hold on
plot(xf,vb1,'Displayname','k_0 = 1')
plot(xf,vb2,'Displayname','k_1 = 12')
plot(xf,vb3,'Displayname','k_3 = 30')
plot(xf,vb4,'Displayname','k_4 = 40')
plot(xf,vb5,'Displayname','k_5 = 85')
plot(xf,vb6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('v_b','Fontsize',16)
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;

%% delta
figure
hold on
plot(xf,delta1,'Displayname','k_0 = 1')
plot(xf,delta2,'Displayname','k_1 = 12')
plot(xf,delta3,'Displayname','k_3 = 30')
plot(xf,delta4,'Displayname','k_4 = 40')
plot(xf,delta5,'Displayname','k_5 = 85')
plot(xf,delta6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\delta','Fontsize',16)
set(gca, 'YScale', 'log')                       % log scale y-axis
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;

%% delta_b
figure
hold on
plot(xf,deltab1,'Displayname','k_0 = 1')
plot(xf,deltab2,'Displayname','k_1 = 12')
plot(xf,deltab3,'Displayname','k_3 = 30')
plot(xf,deltab4,'Displayname','k_4 = 40')
plot(xf,deltab5,'Displayname','k_5 = 85')
plot(xf,deltab6,'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\delta_b','Fontsize',16)
set(gca, 'YScale', 'log')
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;

%% abs(delta_b)
figure
hold on
plot(xf,abs(deltab1),'Displayname','k_0 = 1')
plot(xf,abs(deltab2),'Displayname','k_1 = 12')
plot(xf,abs(deltab3),'Displayname','k_3 = 30')
plot(xf,abs(deltab4),'Displayname','k_4 = 40')
plot(xf,abs(deltab5),'Displayname','k_5 = 85')
plot(xf,abs(deltab6),'Displayname','k_6 = 100')
xlabel('x','Fontsize',16)
ylabel('\delta_b','Fontsize',16)
set(gca, 'YScale', 'log')
lgd = legend('show','Location','Northwest');
lgd.FontSize = 12;






