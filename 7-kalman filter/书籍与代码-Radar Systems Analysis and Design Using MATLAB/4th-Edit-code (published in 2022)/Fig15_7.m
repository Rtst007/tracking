% Use this program to reproduce Fig. 15.7 of text
clc
clear all
close all
eps =0.00;
N = 8;
rect(1:32) = 1;
ham = hamming(32);
han = hanning(32);
blk = blackman(32);
k3 = kaiser(32,3);
k6 = kaiser(32,6);
RECT = 20*log10(abs(fftshift(fft(rect, 1024)))./32 +eps);
HAM =  20*log10(abs(fftshift(fft(ham, 1024)))./32 +eps);
HAN =  20*log10(abs(fftshift(fft(han, 1024)))./32+eps);
BLK = 20*log10(abs(fftshift(fft(blk, 1024)))./32+eps);
K6 = 20*log10(abs(fftshift(fft(k6, 1024)))./32+eps);
x = linspace(-1,1,1024);
figure
plot(x,RECT,'k:',x,HAM,'k--',x,HAN,'k','linewidth',1.);
xlabel('x')
ylabel('Window')
grid
axis tight
legend('Rectangular','Hamming','Hanning')
figure
plot(x,RECT,'k:',x,BLK,'k--',x,K6,'k','linewidth',1.)
xlabel('x')
ylabel('Window')
legend('Rectangular','Blackman','Kasier at \beta = 6')
grid
axis tight

