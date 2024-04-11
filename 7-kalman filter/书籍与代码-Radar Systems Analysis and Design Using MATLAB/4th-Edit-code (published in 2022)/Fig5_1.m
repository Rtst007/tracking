clc
clear all 
close all 
tau= 2;
B = 5;
fs = 2.5;
range= [200 500];
Rwin =750;
k_noise = 0;
n_tgt =  length(range); 
tgt_rcs = [ 1 1.];
c = 3.0e8;
tmax = 2*Rwin/c; 
tau= tau*1e-6; 
B = B *1e6;
fs = fs* 2*B;
ts=	1/fs; %sampling interval
n =	ceil(tau*fs) ; % # of smaples 
t = 0:ts:ts*(n-1);
x = exp(-i*pi*(B/tau) .* t.^2);
xn = x +k_noise*.5*randn(1,n);
t_delay= 2*range/c;	%time delay for target within reecieve window y = zeros(n_tgt,ceil(max(t_delay)/ts)+n-1);
for nt = 1:n_tgt
y(nt,ceil(max(t_delay(nt)/ts)):ceil(max(t_delay(nt)/ts))+n-1) = ...
    tgt_rcs(nt) * real(xn);
end
t_rec = 0:ts:ts*(ceil(max(t_delay)/ts) + n-2); 
xn_ttl = ones(1,n_tgt) * y;
% MF Computation
h = ifft(conj(fft(real(x)))); % replica
% signal processing of all targets
% perform pulse compression 
y = conv(h,xn_ttl);
t_out = 0:ts:ts*(length(y)-1);
t_out = t_out- 2e-6; 
length(t_out);
figure 
t_out = t_out + 2e-6; 
r=	(c/2) .*(t_out - tau)+2.5;
% subplot(2,1,2) 
plot(r,20*log10(abs(y)+.0001),'b','linewidth', 1.5);
grid on
xlim([0 Rwin]);
xlabel('Target range within receive window in meters') 
ylim([min(y)-25 max(y)+25]);
ylabel ('Receiver output vs. range in dB') 
%title(['2 equal RCS targets within receive window @ ',num2str(range),' meters']) 
title(['2 targets within a receive window @ ',num2str(range),' meters'])
%end
