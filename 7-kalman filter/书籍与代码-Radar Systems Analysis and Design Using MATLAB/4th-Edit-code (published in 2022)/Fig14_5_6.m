% Generates plot like Fig. 14.5 and Fig. 14.6 
% Enter scatterer spacing, in meters
clc
close all
clear all
distance = input('Enter scatterer spacing, in meters \n');
% Enter frequency band
freqli = 8;%input('Enter lower frequency limit in Hz \n');
frequi = 12;%input('Enter upper frequency limit in Hz \n');
freql = freqli*1e9;
frequ = frequi*1e9;
[rcs] = rcs_frequency (distance, frequ, freql);
N = size(rcs,2) ;
freq = linspace(freql,frequ,N)./1e9;
figure (1);
plot(freq,rcs,'b', 'linewidth',1);

grid on;
xlabel('\bfFrequency');
ylabel('\bfRCS in dBsm');
title({['\bf frequency =  ',num2str(freqli),' to ', num2str(frequi),' GHz','\bf , spacing = ',num2str(distance), ' m']})

