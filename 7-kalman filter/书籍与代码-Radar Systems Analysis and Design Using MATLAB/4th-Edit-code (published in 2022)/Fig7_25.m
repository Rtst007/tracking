%
clc
close all
clear all
f0 = 16e9;

v = -5000;
b =100e6;
range = 300,
tau=2.5e-6,
nnoise =0.1
rcs = 1;
fs =5*b; 
n = ceil(tau*fs); %Number of sample points of chirp signal
c = 3e8;
lambda = c/f0;
ts = 1/fs;
t = 0:ts:ts*(n-1);
figure
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
tout = 0:ts:ts*(length(filter_out)-1);

delr = (tout)*c/2;

 
H =  conj(fft(replica,nfft));
h= (b*tau/n)*ifft((H.*fft(ycomp)));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
filter_out = filter_out./max(filter_out);
subplot(2,1,1)
plot(delr,20*log10( filter_out+0.00004));
xlabel('\bf Target range within receive window in meters');
ylabel('\bf Normalized MF output in dB')
grid

%*************************************

x1 =  exp(i * pi * (b/tau) .* t.^2);
replica = exp(i * pi * (b/tau) .* t.^2);
ntarget = size(range,2);
time_delay = 2*(range+(tau*v*f0/b))/c;
y1 = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
for k = 1:ntarget
    y1(k,ceil(time_delay(k)/ts):ceil(time_delay(k)/ts)+n-1) = x1*rcs(k);
end
ycomp = ones(1,ntarget)*y1;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));

tsamp = 0:ts:ts*ceil((length(ycomp)-1));

nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);

delr = (tout)*c/2;
 
H =  conj(fft(replica,nfft));
h= (b*tau/n)*ifft((H.*fft(ycomp)));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
filter_out = filter_out./max(filter_out);
subplot(2,1,2)
plot(delr,20*log10( filter_out+0.00004));
xlabel('\bf Target range within receive window in meters');
ylabel('\bf Normalized MF output in dB')
grid





