function [resp] = single_canceler (fofr1)
% single delay canceller
eps = 0.00001;
fofr = 0:0.01:fofr1;
arg1 = pi .* fofr;
resp = 4.0 .*((sin(arg1)).^2);
max1 = max(resp);
resp = resp ./ max1;
figure('units', 'normalized', 'position', [0.1, 0.10, 0.3, 0.4]);
subplot(2,1,1)
plot(fofr,resp,'k','linewidth',1)
xlabel ('\bf Normalized frequency in f/fr')
ylabel( '\bf Amplitude response in Volts')
grid
subplot(2,1,2)
resp=10.*log10(resp+eps);
plot(fofr,resp,'k','linewidth',1);
axis tight
grid
xlabel ('\bf Normalized frequency in f/fr')
ylabel( '\bf Amplitude response in dB')
end

