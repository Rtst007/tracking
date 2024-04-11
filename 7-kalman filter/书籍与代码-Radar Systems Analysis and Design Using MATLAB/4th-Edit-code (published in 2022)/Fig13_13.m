% Use this program to reproduce Fig. 13.13 of text
clc
clear all
close all
index = 0.;
for pd = 0.01:.025:1
    index = index + 1;
    [Lf,Pd_Sw5] = fluct_loss(pd, 1e-9,1,1); 
    Lf1(index) = Lf;
    [Lf,Pd_Sw5] = fluct_loss(pd, 1e-9,1,4);
    Lf4(index) = Lf;
    
end
pd = 0.01:.025:1;
figure (3)
semilogy(Lf1,pd,'k',Lf4, pd,'K-.','linewidth',1.5)
ylabel('\bfProbability of detection')
xlabel('\bfFluctuation loss - dB')
legend('Swerling I & II','Swerling III & IV')
title('P_f_a = 10^-^9, n_p = 1')
grid on




%[Lf,Pd_Sw5] = fluct_loss(pd, pfa, np, sw_case)