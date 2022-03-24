% Program to make Moire patterns with cos lines.
clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%% Creating a random Moire Pattern
% I0=randn,  A=randn+1i*randn                     % generate random parameters
% wx0=rand*2*pi-pi; wy0=rand*2*pi-pi;
% %wx=0.880134906842277;wy=0.880134906842277;     % value for 45 degree
% H=256;K=512;                                    % size of image
% x=1:H; x=repmat(x',[1,K]); 
% y=1:K; y=repmat(y,[H,1]);
% moireImage=I0+real(A*exp(1i*wx0*x+1i*wy0*y));   % create the image
% ft = fftshift(fft2(moireImage));                % Fourier Transform of Moire Pattern

%% Load an image
moireImage0= imread("moire1.jpg");            
imshow(moireImage0); 
[D,E,~]=size(moireImage0);

%% j
H=10;K=10; moireImage=zeros(H,K);
moireImage0=rgb2gray(moireImage0);
moireImage=imcrop(moireImage0,[0 0 H K]);
imshow(moireImage)
ft=fftshift(fft2(moireImage));                       % Fourier Transform of Moire Pattern
abs(ft)/1000;
%imshow(abs(ft)/100000)
%% Plotting the original image and FT of it
    
figure()
subplot(2, 2, 1); imshow(moireImage)           % Moire Pattern

subplot(2,2,2); imshow(abs(ft))                  % fft2 with fftshift

subplot(2,2,3); imshow(abs(ft), [0, 300])        % Choosing a suitable value to show the maxima

subplot(2,2,4); imshow(abs(ft), [0 10000])

%% Finding the location of the maximum point

J = abs(ft)/1000;               % the adjusted image
if mod(H,2)==0
    xorigin = H/2;                 % locate the origin
else
    xorigin = H/2+0.5;
end
if mod(K,2)==0
    yorigin = K/2
else
    yorigin=K/2+0.5;
end
J(xorigin,yorigin) = 0;           % Set the magnitude of origin to 0

[J0, jx] = max(J);              % locate the maxima
[J1, jy] = max(pagetranspose(J));

[~,jx] = max(J1);
[~,jy] = max(J0);

wx = xorigin - jy;
wy = yorigin - jx;

wx = 2 * pi * wx/H;              % calculate wx and wy
wy = 2 * pi * wy/H;

%% Computing X, the matrix containing I0 and A
A = ones(H * H,1);                          % creating empty matrixs
B = zeros(H);
C = zeros(H);
for x=1:H
    for y=1:K
        pow = 1i * wx * x + 1i * wy * y; 
        B(x,y) = exp(pow);
        C(x,y) = exp(-pow); 
    end
end
I = double(reshape(moireImage,[],1));
I=I./max(I);
M = [A,B(:),C(:)];                          % Put the vector in one matrix
X = M\I(:)

% %% Regenerate the original Moire pattern with matrix X
% I2=X(1,1);  A2=2*X(2,1);                             % Get the parameters from matrix X
% x=1:H; x=repmat(x',[1,K]); 
% y=1:K; y=repmat(y,[H,1]);
% moireImage2=I2+real(A2*exp(1i*wx*x+1i*wy*y));     % create the image
% ft2 = fftshift(fft2(moireImage2));                  % Fourier Transform of Moire Pattern
% 
% %% Plotting the regenerated image and FT of it
%     
% figure();
% subplot(2, 2, 1); imshow(moireImage2);           % Moire Pattern
% 
% subplot(2,2,2); imshow(abs(ft))                  % fft2 with fftshift
% 
% subplot(2,2,3); imshow(abs(ft), [0, 300])        % Choosing a suitable value to show the maxima
% 
% subplot(2,2,4); imshow(abs(ft), [0 10000])