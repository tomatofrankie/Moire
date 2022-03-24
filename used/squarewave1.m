clc
clear all
close all

x = linspace(-4*pi, 4*pi, 500);
y = 1/2 + 2/pi * sin(x) + 2/(3*pi) * sin(3*x)+ 2/(5*pi)*sin(5*x);
axis on
plot(x, y, "Color",[0,0,0])
