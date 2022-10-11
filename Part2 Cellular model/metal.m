%% MT
OD = 1;
ODunit = 8E+8;
b1 = 0.084;% the maximum speed of hill equation uM*min-1   
Kx1 = 0.00005;%uM
Kr1 = 1;% uM
LacItot = 0.01;% may be related to the numebr of E.coil  uM the other data is 498?
IPTG = 1000*0.001/OD/ODunit*0.92/60/1E-15*1E-14/(0.001/OD/ODunit);% induce
n = 2;% hill cofficient
rm = 3.4E-3 * 60;%min-1
rp = 4.88E-5 * 60;%min-1
ktr = 0.44 * 60;%min-1
kt = 1.8E-4 * 60; % min-1
[t,y] = ode45(@(t,y) MBPsyn3_2(b1,Kx1,Kr1,LacItot,IPTG,n,rm,rp,ktr,y,kt),0:4320,[0,0,0]);
b = y(:,3)*1E-15;

plot(t,b,"LineWidth",1.1);
hold on;

%%  MBP
b1 = 0.084;% the maximum speed of hill equation uM*min-1   
Kx1 = 0.00005;%uM
Kr1 = 1;% uM
LacItot = 0.01;% may be related to the numebr of E.coil  uM the other data is 498?
IPTG = IPTG*1000*0.001/OD/ODunit*0.92/60/1E-15*1E-14/(0.001/OD/ODunit);% induce
n = 2;% hill cofficient
rm = 4.4E-3 * 60;%min-1
rp = 6.3E-5 * 60;%min-1
ktr = 0.57 * 60;%min-1
kt = 1.8E-4 * 60; % min-1
[t,y] = ode45(@(t,y) MBPsyn3_2(b1,Kx1,Kr1,LacItot,IPTG,n,rm,rp,ktr,y,kt),0:4320,[0,0,0]);
c = y(:,3)*1E-15;

plot(t,c,"LineWidth",1.1);
title("The expression of CsgA·MBP/MT in one cell","FontSize",16);
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
ylabel("Amount of substance(μmol)");
xlabel("time(min)");
legend("CsgA·MT","CsgA·MBP");
grid on;