OD = 0.7;
ODunit = 8E+8;% the number of cells/ml

Kx1 = 0.00005;% unit:μM
Kr1 = 1;% unit:μM
LacItot = 0.01;% the normal concentration of LacI in one E.coli,unit:μM
n = 2;% hill cofficient
IPTG = (0:0.05:2);% the range of IPTG,unit:mM
IPTG1 = IPTG*1000*0.001/OD/ODunit*0.92/60/1E-15*1E-14/(0.001/OD/ODunit);% The concentration of IPTG in the cell
c=0.0014./ (1 + ((LacItot./Kx1) .*(1 ./(1 + (IPTG1.^n)./Kr1))));% promoter activity, unit:μM/s

plot(IPTG,c,'ok-', 'linewidth', 1.1, 'markerfacecolor', "c",'MarkerSize',8);
title("The relationship between the concentration of IPTG and promoter activity","FontSize",16,'Color',	'#D95319');
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
xlabel("The concentration of IPTG(mM)","FontSize",16,'Color','#D95319');
ylabel("Promoter activity(μM/s)","FontSize",16,'Color',	'#D95319');
legend("tac promoter activity");
grid on;
box on;