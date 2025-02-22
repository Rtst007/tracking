% Use this program to reproduce Fig. 3.5 of text.
close all
clear all
pt = 1.5e+6; % peak power in Watts
freq = 5.6e+9; % radar operating frequency in Hz
g = 45.0; % antenna gain in dB
sigma = 0.1; % radar cross section in m squared
b = 5.0e+6; % radar operating bandwidth in Hz
nf = 3.0; %noise figure in dB
loss = 6.0; % radar losses in dB
np = 1;
range = linspace(25e3,225e3,500); % range to target from 5 Km 225 Km, 1000 points
snr1 = lprf_req(pt, g, freq, sigma, np, b, nf, loss, range);
snr2 = lprf_req(pt, g, freq, sigma, 5*np, b, nf, loss, range);
snr3 = lprf_req(pt, g, freq, sigma, 20*np, b, nf, loss, range);
% plot SNR versus range
figure(1)
subplot(2,1,1)
rangekm  = range ./ 1000;
semilogx(rangekm,snr3,'k',rangekm,snr2,'k:',rangekm,snr1,'k-. ','linewidth',1.5)
grid
legend('n_p = 20','n_p = 5','n_p = 1')
xlabel ('Detection range - Km');
ylabel ('SNR - dB');
np = linspace(1,500,500);
range = 150e3;
snr1 = lprf_req(pt, g, freq, sigma, np, b, nf, loss, range);
snr2 = lprf_req(pt, g, freq, 10*sigma, np, b, nf, loss, range);
subplot(2,1,2)
semilogx(np,snr2,'k',np,snr1,'k: ','linewidth',1.5)
grid
legend('\sigma = 0 dBsm','\sigma = -10 dBsm')
xlabel ('No. of pulses');
ylabel ('SNR - dB');




