clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%%
x = linspace(-4*pi, 4*pi, 500);
subplot(1,2,1)
plot(x,1/2 + square(x)/2, "Color",[0,0,0])
%%
ft = fftshift(fft(x));
%subplot(1,2,2)
plot(x-0.025,abs(ft), "Color", [0,0,0])