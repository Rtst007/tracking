function [resp] = double_canceler(fofr1)
eps = 0.00001;
fofr = 0:0.01:fofr1;
arg1 = pi .* fofr;
resp = 4.0 .* ((sin(arg1)).^2);
max1 = max(resp);
resp = resp ./ max1;
resp2 = resp .* resp;
figure('units', 'normalized', 'position', [0.1, 0.10, 0.3, 0.4]);
subplot(2,1,1);
plot(fofr,resp,'k--',fofr, resp2,'k','linewidth',1);
xlabel ('\bf Normalized frequency f/fr')
ylabel ('\bf Amplitude response - Volts')
resp2 = 20. .* log10(resp2+eps);
resp1 = 20. .* log10(resp+eps);
grid on
subplot(2,1,2)
plot(fofr,resp1,'k--',fofr,resp2,'k');
legend ('\bf Single canceler','double canceler')
xlabel ('\bf Normalized frequency f/fr')
ylabel ('\bf Amplitude response - dB')
grid on
end