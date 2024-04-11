function [y] = stretch(nscat,taup,f0,b,rrec,scat_range,scat_rcs,winid)
% This function implements the stretch processor
%
% Inputs
    % nscat         == number of point scatterers within the receive window	
    % taup          == uncompressed pulse width	in seconds
    % f0            == chirp start frequency in Hz
    % b             == chirp bandwidth in Hz
    % scat_range	== vector of scatterers’ range in m
    % rrec          == range receive window in m
    % scat_rcs      == vector of scatterers’ RCS in m^2
    % win           == 0 = no window; 1 = Hamming; 2 = Kaiser with parameter pi; 3 = Chebychev side-lobes at -60dB	
% Outputs
    % y             == compressed output in volts
%
eps = 1.0e-14;
htau = taup / 2.;
c = 3.e8;
trec = 2. * rrec / c;
n = fix(2. * trec * b);
m = power_integer_2(n);
nfft = 2.^m;
x(nscat,1:n) = 0.;
y(1:n) = 0.;
if( winid == 0.)
   win(1:n) = 1.;
   win =win';
else
   if(winid == 1.)
      win = hamming(n);
   else
      if( winid == 2.)
         win = kaiser(n,pi);
      else
         if(winid == 3.)
            win = chebwin(n,60);
         end
      end
   end
end
deltar = c / 2. / b;
max_rrec = deltar * nfft / 2.;
maxr = max(scat_range);
if(rrec > max_rrec | maxr >= rrec )
   'Error. Receive window is too large; or scatterers fall outside window'
   return
end
t = linspace(0,taup,n);
for j = 1:1:nscat
    range = scat_range(j);% + rmin;
   psi1 = 4. * pi * range * f0 / c - ...
      4. * pi * b * range * range / c / c/ taup;
   psi2 = (4. * pi * b * range / c / taup) .* t;
   x(j,:) = scat_rcs(j) .* exp(2*i * psi1 + 2*i .* psi2);
   y = y + x(j,:);
end
y =y +0.25*randn(1,size(y,2));
figure(1)
time = linspace(0,trec,n);
plot(time*1e6,20*log10(abs(y)),'b','linewidth',.75)
xlabel ('\bf Relative delay within T_r_e_c - \mu s')
ylabel ('\bf Noisy uncompressed echo')
grid
ywin = y .* win';
yfft = fft(y,n) ;
out= (abs(yfft));
figure(2)
delinc = rrec/ n;
%dist = linspace(-delinc-rrec/2,rrec/2,n);
dist = linspace(0, rrec,n);
out = out ./ max(abs(out));
plot(dist,20*log10(out+eps),'b','linewidth',1)
xlabel ('\bf Relative range in meters')
ylabel ('\bf Normalized pulsed copressed signal in dB')
axis auto
grid

