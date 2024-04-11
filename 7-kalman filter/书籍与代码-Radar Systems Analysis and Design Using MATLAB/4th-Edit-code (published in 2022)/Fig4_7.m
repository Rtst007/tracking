% this program reproduces Fig. 4.7 of text
clc
close all
clear all
Rmax = 1600;        % Km
el = [1 2 5 10];   % elevation angle in deg
H = 78.11;          % km
No = 313;           % refractivity at earth surface
Ce = 0.1439;        % km^-1
pmax = 1.25e5;      % maximum electron density at hm
hm = 300.73;        % % height for maximum electron contents in Kmkm
f = 450e6;          % hz, center frequency
[deltaR, Rm, Rt] = refraction(Rmax, el(1),H, No, Ce, pmax, hm, f);
figure
plot(cumsum(Rm), deltaR .*1000, 'b--','linewidth',1.5)
hold on
[deltaR, Rm, Rt] = refraction(Rmax, el(2),H, No, Ce, pmax, hm, f);
plot(cumsum(Rm), deltaR .*1000, 'b:','linewidth',1.5)
hold on
[deltaR, Rm, Rt] = refraction(Rmax, el(3),H, No, Ce, pmax, hm, f);
plot(cumsum(Rm), deltaR .*1000, 'b -.','linewidth',1.5)
hold on
[deltaR, Rm, Rt] = refraction(Rmax, el(4),H, No, Ce, pmax, hm, f);
plot(cumsum(Rm), deltaR .*1000, 'b','linewidth',1.5)
hold off
grid on
xlabel('\bfRange - Km ')
ylabel('\bfRange Error - meters')
legend('\beta = 1 deg','\beta = 2 deg','\beta = 5 deg','\beta = 10 deg')
title('Frequency f = 9.5 GHz')