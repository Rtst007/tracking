% generates Fig. 14.3 of text
clc
close all
clear all
% Enter scatterer spacing, in meters
distance = .2; %input('Enter scatterer spacing, in meters \n');
% Enter frequency
freqi = 5.6;%input('Enter Enter frequency in Hz \n');
freq = freqi *1e9;
rcs = rcs_aspect(distance,freq);
figure (1);
aspect_degrees = 0.:.05:180.;
plot(aspect_degrees,rcs,'b','linewidth',1);
grid;
xlabel('\bfaspect angle - degrees');
ylabel('\bfRCS in dBsm');
title({['\bf frequency =  ',num2str(freqi),' GHz','\bf , spacing = ',num2str(distance), 'm']})

