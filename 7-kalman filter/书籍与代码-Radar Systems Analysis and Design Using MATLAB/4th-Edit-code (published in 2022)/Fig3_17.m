% This program generates Fig. 3.17 of text
clc
clear all
close all
pt =10.0e+3; pt_db = 10*log10(pt);
g = 35.0;       % antenna gain in dB
freq = 5.6e+9;  lambda = 3e8 / freq;
lambda_db = 10*log10(lambda^2);
sigma = 10 ;  
b = 500.0e+3;  b_db = 10*log10(b);
range = 22.2  
range_db = 10*log10(range * 1000.);
gprime = 10.0;   sigmadb = 10*log10(sigma);
loss = 0.01;    
rangej = 22.2; rangej_db = 10*log10(rangej * 1000.);
pj = 1000;   pj_db = 10*log10(pj);  
bj = b;%10.0e+6; 
bj_db = 10*log10(bj);
gj = 30.0; 	
lossj =0.3;
factor = 10*log10(4.0 *pi);
%[BR_range] = soj_req (pt, g, sigma, b, freq, loss, range,pj, bj,gj, lossj, gprime, rangej)*1000

BR_range = ((pt * 10^(2.0*g/10) * sigma * bj * 10^(lossj/10) * (rangej*1000)^2) / (4.0 * pi * pj * 10^(gj/10) * 10^(gprime/10) * ...
   b * 10^(loss/10)))^.25 

% soj_req (pt, g, sigma, b, freq, loss, range, pj, bj,gj, lossj, gprime, rangej);
 s_at_br = pt_db + 2.0 * g + lambda_db + sigmadb - 3.0 * factor - 4.0 * 10*log10(BR_range) - loss  
   soj1 = s_at_br;%bj_db+gj+gprime+lambda_db+b_db-lossj-2*factor-2*rangej_db-bj_db-loss
index =0;
for ran_var = .01:.5:1000;
   index = index + 1;
   ran_db = 10*log10(ran_var * 1000.0);
  s(index) = pt_db + 2.0 * g + lambda_db + sigmadb - 3.0 * factor - 4.0 * ran_db - loss;% + s_at_br;
  soj(index) = soj1;%bj_db+gj+gprime+lambda_db+b_db-lossj-2*factor-2*rangej_db-bj_db-loss;
   % soj(index) = s_at_br ;%- s_at_br;
end

ranvar = .01:.5:1000;
%ranvar = ranvar ./BR_range;
semilogx (ranvar,s,'b',ranvar,soj,'r','linewidth',1.5);
axis([.1 1000 -150 0])
xlabel ('Range normalized to cross-over range');
legend('Target echo','SOJ')
ylabel ('Relative signal or jamming amplitude - dB');
grid
