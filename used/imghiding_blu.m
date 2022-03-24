clc;
close all;  
clear;  
s1=512;s2=512;                                  % size of the image
x=repmat(((1:s2)-(s2+1)/2)/(s2-1),[s1,1]);      % x,y variables uniform in
y=-repmat(((1:s1)'-(s1+1)/2)/(s1-1),[1,s2]);    % the interval [-0.5,0.5]
 
M=@(p,lambda)(1+cos(2*pi./lambda.*p))/2;
lambda=@(y)0.1*(y>=0)+0.18*(y<0);               % image is made of two different frequencies (0.1 and 0.18)
a0=fzero(@(x)besselj(0,x),2);                   % choice of a zero of J_0(x)
a=@(x)a0/(2*pi/0.18);                           % targeting one of the two frequencies (0.18) to be cancelled
w=2*pi;phi=0;u=@(x,t)a(x).*sin(w*t+phi);
 
K=100;I=0*x;
for k=1:K
    t=k/K;
    I=I+M(angle(x+i*y)-u(x,t),lambda(y));
end
I=I/100;
imshow(I,[min(I(:)),max(I(:))])