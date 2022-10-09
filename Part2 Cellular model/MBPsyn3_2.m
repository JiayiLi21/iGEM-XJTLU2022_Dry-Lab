function dydt = MBPsyn3_2(b1,Kx1,Kr1,LacItot,IPTG,n,rm,rp,ktr,y,kt)
dydt = zeros(3,1);
dydt(1) = b1 ./ (1 + ((LacItot./Kx1) .*(1 ./(1 + (IPTG.^n)./Kr1)))) - rm*y(1);
dydt(2) = ktr*y(1) - rp*y(2) - kt*y(2);
dydt(3) = kt*y(2);
end