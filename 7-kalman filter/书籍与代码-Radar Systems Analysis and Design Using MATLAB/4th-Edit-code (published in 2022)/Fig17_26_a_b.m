
%Generates Fig 17.26-a and 17.26-b
clear all
eps = 0.0000001;
npts = 5000;
del = 1./ 5000.;
t = 0. : del : 1.;
% generate input sequence
inp = 1.+ t.^3 + .5 .*t.^2 + cos(2.*pi*5 .* t) ;
% read the intial estimate for the state vector
X0 = [2,.1,.01]';
% this is the update interval in seconds
T = 100. * del;
% this is the value of the smoothing coefficient
smoocof = .1;
nvar = eps;
[residual, estimate] = ghk_tracker (X0, smoocof, inp, npts, T, nvar);

figure(1)
NN = 4999.;
n = 1:NN;
subplot(2,1,1)
plot(n,estimate(1:NN),'k')
xlabel ('Sample number')
ylabel ('Predicted position')
grid on

subplot(2,1,2)
plot(n,inp(1:NN),'k')
xlabel ('Sample number')
ylabel ('Truth position')
grid on


figure(3)
plot (residual(1:150),'k')
xlabel ('Sample number')
ylabel ('Residual error')
grid
