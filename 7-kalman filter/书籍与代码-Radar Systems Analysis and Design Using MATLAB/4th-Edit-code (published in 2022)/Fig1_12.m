% This function implements the matched filter processor
%
% Inputs
    % tau          == uncompressed pulse width in seconds
    % b            == chirp bandwidth in Hz
    % range        == scatterers’ relative range in m
    % rcs          == vector of scatterers’ RCS in meter squared
% Output
    % y             == normalized compressed output
%
clc
close all
clear all
tau = 1e-6;
b = 50e6;
fs =5*b; 
rcs = 1;
range = 150;
nnoise = 0.05;
n = ceil(tau*fs); %Number of sample points of chirp signal
c = 3e8;
ts = 1/fs;
t = 0:ts:ts*(n-1);

% Transmitted LFM Chirp pulse:
x =  exp(i * pi * (b/tau) .* t.^2);
replica = exp(i * pi * (b/tau) .* t.^2);

ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
for k = 1:ntarget
    y(k,ceil(time_delay(k)/ts):ceil(time_delay(k)/ts)+n-1) = x*rcs(k);
end
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));

tsamp = 0:ts:ts*ceil((length(ycomp)-1));



nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0.0:ts:ts*(length(filter_out)-1);
% subplot(2,2,4); 
% plot(tout*1e6,20*log10( filter_out+0.0004));
% title('\bf Output of Matched Filter');
% xlabel('\bf Time \mu s');
% grid
% ylabel ('\bf Pulse compression gain in dB')

delr = (tout)*c/2;
figure
subplot(2,1,1)
filter_out = filter_out./max(filter_out);
plot(70*tout*1e6,20*log10( filter_out+0.0004),'linewidth',1.5);
title('\bf Output of Matched Filter');
xlabel('\bf Time \mu s');
grid
ylabel ('\bf Output vs. time')
title({['\bf Targets range 10.5 km', '; delay = 70 \mu sec']})

% determine Hamming window

H =  conj(fft(replica,nfft));
h= (b*tau/n)*ifft((H.*fft(ycomp)));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
filter_out = filter_out./max(filter_out);
subplot(2,1,2)
plot((10350+delr)./1000,20*log10( filter_out+0.00004));
xlabel('\bf Target range in km');
ylabel('\bf Output vs. range')
grid


