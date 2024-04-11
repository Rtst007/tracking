% Use this program to reproduce Fig. 8.20 of the textbook.
clc
clear all
close all
taup = 0.4;
pri = 1;
n = 3;
bw = 50;

xlfm = train_ambg_lfm(taup, n, pri, bw);
doppler = linspace(-1/taup,1/taup,size(xlfm,1));

timelfm = linspace(-taup,taup,size(xlfm,2));

figure(1)
contour(timelfm, doppler, (xlfm));
%surf(time, doppler, x); shading interp; view(0,90);
xlabel('\bf Delay in seconds');
ylabel('\bf Doppler in Hz');
grid;
axis tight;
title('\bf LFM pulse train, B\tau = 40, N = 3 pulses')

xpt = train_ambg(taup, n, pri);
timept = linspace(-n*taup,n*taup,size(xpt,2));
doppler = linspace(-1/taup,1/taup,size(xpt,1));

figure(2)
contour(timept, doppler, (xpt));
xlabel('\bf Delay in seconds');
ylabel('\bf Doppler in Hz');
grid;
axis tight;
