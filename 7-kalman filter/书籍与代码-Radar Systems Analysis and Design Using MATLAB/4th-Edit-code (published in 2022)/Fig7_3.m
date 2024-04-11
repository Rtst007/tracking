% use this program to reproduce Fig. 7.3 of text
clear all
close all
clc
nscat = 2; %two point scatterers
tau = 10e-6; % 10 micro second pulse
b = 10.0e6; % 50 MHz bandwidth
range = [1000 1500] ; % scatterers are 15 and 25 meters into window
rcs = [1 2]; % RCS 1 m^2 and 2m^2
nnoise = 01; %no window used
fs =2.5*b; 
n = ceil(tau*fs); %Number of sample points of chirp signal
c = 3e8;
ts = 1/fs;
t = 0:ts:ts*(n-1);

% Transmitted LFM Chirp pulse:
x =  exp(i * pi * (b/tau) .* t.^2);
replica = exp(i * pi * (b/tau) .* t.^2);
figure
subplot(2,1,1)
% plot(t*1e6,real(x));
% xlabel ('\bf time in \mu s')
% ylabel ('\bf Real part of replica')
% subplot(2,2,2)
% plot(t*1e6,imag(x));
% xlabel ('\bf time in \mu s')
% ylabel ('\bf Imaginary part of replica')
ntarget = size(range,2);
time_delay = 2*range/c;
y = zeros(ntarget,ceil(max(time_delay/ts))+n-1);
for k = 1:ntarget
    y(k,ceil(time_delay(k)/ts):ceil(time_delay(k)/ts)+n-1) = x*rcs(k);
end
ycomp = ones(1,ntarget)*y;
ycomp =ycomp +nnoise*randn(1,size(ycomp,2));

tsamp = 0:ts:ts*ceil((length(ycomp)-1));


plot(tsamp*1e6,ycomp,'linewidth',1);
grid; 
title({['\bf Targets range [',num2str(range),'] meters'];...
    ['\bf Traget RCS  [',num2str(rcs), '] m^2']})
xlabel('\bf Time in \mu s');
ylabel ('\bf Uncompressed signal ')


nfft =size(ycomp,2)
H =  (conj(fft(replica,nfft)));
h = (b*tau/n)*ifft(H.*fft(ycomp));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
subplot(2,1,2); 
% 
 delr = (tout)*c/2+(c/4/b);
filter_out = filter_out./max(filter_out);

% determine Hamming window
win = hamming(size(replica,2));
 
H =  conj(fft(replica.*win',nfft));
h= (b*tau/n)*ifft((H.*fft(ycomp)));
filter_out = h;
tout = 0:ts:ts*(length(filter_out)-1);
filter_out = (filter_out./b ./tau);
subplot(2,1,2)
plot(delr,10*log10( filter_out+0.00004),'linewidth',1);
xlabel('\bf Target range within receive window in meters');
ylabel('\bf Pulse compressed signal in dB')
grid
axis tight

