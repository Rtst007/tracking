% Generates Fig.s 17.22 
close all
clc
clear all
eps = 0.0000001;
npts = 5000;
del = 1./ 5000.;
t = 0. : del : 1.;
% generate input sequence
inp = 1.+ t.^3 + .5 .*t.^2 ;%+ cos(2.*pi*5 .* t) ;
% read the intial estimate for the state vector
X0 = [2,.1,.01]';
% this is the update interval in seconds
T = 100. * del;
% this is the value of the smoothing coefficient
xi = .1;
[residual, estimate] = ghk_tracker (X0, xi, inp, npts, T, eps);
figure(1)
plot(inp,'k')
xlabel ('Sample number')
ylabel ('Position')
axis tight
grid on