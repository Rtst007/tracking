% generates Fig 9.10 of text
clc
close all
clear all
fofr=0:0.001:1;
f1=4.*fofr;
f2=5.*fofr;
arg1=pi.*f1;
arg2=pi.*f2;
resp1=abs(sin(arg1));
resp2=abs(sin(arg2));
resp=resp1+resp2;
max1=max(resp);
resp=resp./max1;
figure('units', 'normalized', 'position', [0.1, 0.10, 0.3, 0.6]);

subplot(3,1,1)
plot(fofr,resp1,'linewidth',1);
ylabel('\bfFilter response');xlabel('\bff/fr')
grid on
subplot(3,1,2)
plot(fofr,resp2,'linewidth',1);
ylabel('\bfFilter response');xlabel('\bff/fr')
grid on
subplot(3,1,3)
plot(fofr,resp,'linewidth',1);
ylabel('\bfFilter response')
xlabel('\bff/fr')
grid on

