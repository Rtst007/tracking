close all
clear all

pfa = 1e-6;
TNR = -log(pfa);
SNR =linspace(10,25,200);
snr = 10.^(SNR./10);
 N = 16; % number of cells
 T = N*(exp(-log(pfa)/N)-1);
 %Pd_swr1= exp(-TNR./(1 + SNR));		
 Pd12 = exp(-TNR./(1+snr));
 Pd8 = 1./ (1+(T./(N*(1+snr)))).^N;
 %Pd_3_4 = (1 + (2.*SNR.*TNR./(2+SNR).^2)).* exp((-2.*TNR)./(2+SNR));						
 Pd34 = (1 + (2.*snr.*TNR./(2+snr).^2)).* exp((-2.*TNR)./(2+snr));
figure(1)
 plot(SNR,Pd12,'b',SNR, Pd8,'b:','linewidth', 1.5)
 legend('P_d no CFAR','P_d with CFAR')
 grid on
 xlabel('SNR - in dB','fontweight','bold')
 ylabel('P_d','fontweight','bold')
 xticks([10  11 12 13  14 15 16 17 18 19 20 21 22 23 24 25])
 gtext('CFAR Loss','fontweight', 'bold')
 title('M = 16; P_f_a = 10e^-^6; Swerling I and II target fluctuation','fontweight','bold')
 
 figure(2)
 plot(SNR,Pd34,'b',SNR, Pd8,'b:','linewidth', 1.5)
 legend('P_d no CFAR','P_d with CFAR')
 grid on
 xlabel('SNR - in dB','fontweight','bold')
 ylabel('P_d','fontweight','bold')
 xticks([10  11 12 13  14 15 16 17 18 19 20 21 22 23 24 25])
 gtext('CFAR Loss','fontweight', 'bold') 
 title('M = 16; P_f_a = 10e^-^6; Swerling III and IV target fluctuation','fontweight','bold')
