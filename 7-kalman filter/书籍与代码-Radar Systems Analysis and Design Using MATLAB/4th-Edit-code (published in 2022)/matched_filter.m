function [y] = matched_filter(tau,b,range,rcs,nnoise);
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
fs =2.5*b; 
n = ceil(tau*fs); %Number of sample points of chirp signal
c = 3e8;
ts = 1/fs;
t = 0:ts:ts*(n-1);

% Transmitted LFM Chirp pulse:
x =  exp(i * pi * (b/tau) .* t.^2);
replica = exp(i * pi * (b/tau) .* t.^2);
figure
subplot(2,2,1)
plot(t*1e6,real(x));
xlabel ('\bf time in \mu s')
ylabel ('\bf Real part of replica')
subplot(2,2,2)
plot(t*1e6,imag(x));
xlabel ('\bf time in \mu s')
ylabel ('\bf Imaginary part of replica')
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
for k = 1:ntarget
    y(k,ceil(time_delay(k)/ts):ceil(time_delay(k)/ts)+n-1) = x*rcs(k);
end
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));

tsamp = 0:ts:ts*ceil((length(ycomp)-1));


subplot(2,2,3); plot(tsamp*1e6,ycomp);
title('\bf Return Signal'); 
grid; 
xlabel('\bf Time in \mu s');
ylabel ('\bf Noisy composite retrun')


nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
subplot(2,2,4); 
plot(tout*1e6,20*log10( filter_out+0.0004));
title('\bf Output of Matched Filter');
xlabel('\bf Time \mu s');
grid
ylabel ('\bf Pulse compression gain in dB')

delr = (tout)*c/2+c/4/b;
figure
subplot(2,1,1)
filter_out = filter_out./max(filter_out);
plot(delr,20*log10( filter_out+0.0004));
xlabel('\bf Target range within receive window in meters - No window');
ylabel('\bf Normalized matched filter output in dB')
grid
title({['\bf Targets range [',num2str(range),'] meters'];...
    ['\bf Traget RCS  [',num2str(rcs), '] m^2']})

% determine Hamming window
win = hamming(size(replica,2));
 
H =  conj(fft(replica.*win',nfft));
h= (b*tau/n)*ifft((H.*fft(ycomp)));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
filter_out = filter_out./max(filter_out);
subplot(2,1,2)
plot(delr,20*log10( filter_out+0.00004));
xlabel('\bf Target range within receive window in meters - Window applied');
ylabel('\bf Normalized matched filter output in dB')
grid


