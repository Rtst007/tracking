clc
close all
close all
tau = 1e-6;
b = 25e6;
nnoise =0.1
rcs =1;
range = 150;
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
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
y(1,ceil(time_delay/ts):ceil(time_delay/ts)+n-1) = x*rcs;
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));
tsamp = 0:ts:ts*ceil((length(ycomp)-1));
nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
delr = (tout)*c/2+c/4/b;
plot(delr,20*log10( filter_out+0.00004),'linewidth',1);title('Perfectly matched replica case');
xlabel('range - m');
ylabel('dB')
grid
axis tight
% pulse width mismatch
tau = tau*.9;
x =  exp(i * pi * (b/tau) .* t.^2);
subplot(2,2,2)
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
y(1,ceil(time_delay/ts):ceil(time_delay/ts)+n-1) = x*rcs;
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));
tsamp = 0:ts:ts*ceil((length(ycomp)-1));
nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
delr = (tout)*c/2+c/4/b;
plot(delr,20*log10( filter_out+0.00004),'linewidth',1);
title('10% pulsewdith mismatch case');
xlabel('range - m');
ylabel('dB')
grid
axis tight

% bandwidth mismatch
tau=tau/.9;
b = b*1.1;
x =  exp(i * pi * (b/tau) .* t.^2);
subplot(2,2,3)
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
y(1,ceil(time_delay/ts):ceil(time_delay/ts)+n-1) = x*rcs;
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));
tsamp = 0:ts:ts*ceil((length(ycomp)-1));
nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
delr = (tout)*c/2+c/4/b;
plot(delr,20*log10( filter_out+0.00004),'linewidth',1);
title('10% bandwidth mismatch case');
xlabel('range - m');
ylabel('dB')
grid
axis tight
% pulsewdith & bandwidth mismatch
tau=.9*tau;
x =  exp(i * pi * (b/tau) .* t.^2);
subplot(2,2,4)
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
y(1,ceil(time_delay/ts):ceil(time_delay/ts)+n-1) = x*rcs;
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));
tsamp = 0:ts:ts*ceil((length(ycomp)-1));
nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
delr = (tout)*c/2+c/4/b;
plot(delr,20*log10( filter_out+0.00004),'linewidth',1);
title('10% bandwidth & pulsewidth mismatch case');
xlabel('range - m');
ylabel('dB')
grid
axis tight

