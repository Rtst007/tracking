% Use this program to reproduce Fig. 3.9 of text.
close all
clear all
tsc = 2.5; % Scan time i s2.5 seconds
sigma = 0.1; % radar cross section in m squared
te = 900.0; % effective noise temperature in Kelvins
snr = 15; % desired SNR in dB
nf = 6.0; %noise figure in dB
loss = 7.0; % radar losses in dB
az_angle = 2; % search volume azimuth extent in degrees
el_angle = 2; %serach volume elevation extent in degrees
range = linspace(1e3,1000e3,1000); % range to target from 20 Km 250 Km, 1000 points
pap1 = power_aperture(snr,tsc,sigma/10,range,nf,loss,az_angle,el_angle);
pap2 = power_aperture(snr,tsc,sigma,range,nf,loss,az_angle,el_angle);
pap3 = power_aperture(snr,tsc,sigma*10,range,nf,loss,az_angle,el_angle);
% plot power aperture prodcut versus range
figure(1)
subplot(2,1,1)
rangekm  = range ./ 1000;
semilogx(rangekm,pap1,'b',rangekm,pap2,'b:',rangekm,pap3,'b-.','linewidth',1.5)
grid
legend('\sigma = -20 dBsm','\sigma = -10dBsm','\sigma = 0 dBsm')
xlabel ('Detection range in Km');
ylabel ('Power aperture product in dB');
lambda = 0.03; % wavelength in meters
G = 45; % antenna gain in dB
ae = linspace(1,1000,1000);% aperture size 1 to 25 meter squared, 1000 points
Ae = 10*log10(ae);
range = 250e3; % rnage of interset is 250 Km
pap1 = power_aperture(snr,tsc,sigma/10,range,nf,loss,az_angle,el_angle);
pap2 = power_aperture(snr,tsc,sigma,range,nf,loss,az_angle,el_angle);
pap3 = power_aperture(snr,tsc,sigma*10,range,nf,loss,az_angle,el_angle);
Pav1 = pap1 - Ae;
Pav2 = pap2 - Ae;
Pav3 = pap3 - Ae;
subplot(2,1,2)
semilogx(ae,Pav1,'b',ae,Pav2,'b:',ae,Pav3,'b-.','linewidth',1.5)
grid
xlabel('Aperture size in square meters')
ylabel('Pav in dB')
legend('\sigma = -20 dBsm','\sigma = -10dBsm','\sigma = 0 dBsm')



