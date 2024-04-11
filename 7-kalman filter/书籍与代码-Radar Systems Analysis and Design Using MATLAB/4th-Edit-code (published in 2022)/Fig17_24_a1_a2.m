% Generates Fig.s 17.24a.1 and 17.24-a.2
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
smoocof = .1;
nvar =0.05;
[residual, estimate] = ghk_tracker (X0, smoocof, inp, npts, T, eps);

%ghk_tracker (X0, smoocof, inp, npts, T, nvar);% ghk_tracker (X0, xi, inp, npts, T, eps);

%Figure 17.24-a.1

figure(1)
subplot(2,1,1)
NN = 4999.;
n = 1:NN;
plot(n,estimate(1:NN),'k')
xlabel ('Sample number')
ylabel ('Predicted position')

axis tight
grid on
subplot(2,1,2)
NN = 4999.;
n = 1:NN;
plot(n,inp(1:NN),'k')
xlabel ('Sample number')
ylabel ('Truth psition')
axis tight
grid on

%Figure 17.24-a.b
figure(2)
plot (residual(1:150),'k')
xlabel ('Sample number')
ylabel ('Residual error')
grid on


