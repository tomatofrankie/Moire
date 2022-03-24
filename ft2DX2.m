% Program to get the X matrix of an image
clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%% Load an image
original = rgb2gray(imread("moire1.jpg"));
moireImage0= double(rgb2gray(imread("moire1.jpg")));      % read an image            
[D,E]=size(moireImage0);              % get the size of the image

%% Crop the image into small squares
H=100;K=100;                                  % size of kernel
a=1;b=1;
original = original(a:a+H-1,b:b+K-1);
imshow(original)
moireImage = moireImage0(a:a+H-1,b:b+K-1);        % crop the image to the size of kernel
ft=fftshift(fft2(moireImage));                      % Fourier Transform of kernel

%% Finding the location of the maximum point
J = abs(ft)/1000;                   % the adjusted ft
if mod(H,2)==0                      % locate the origin.
    xorigin = H/2;                  % If the size of H is even
else                                % the origin will be H/2
    xorigin = H/2+0.5;              % else, the origin will be H/2 + 0.5
end
if mod(K,2)==0
    yorigin = K/2;
else
    yorigin=K/2+0.5;
end
J(xorigin+1,yorigin+1) = 0;           % Set the magnitude of origin to 0

[jx, jy] = find(ismember(J, max(J(:))))        % locate the maxima of the ft.

wx = xorigin-jx(1);                            % find the difference between
wy = yorigin-jy(1);                            % the origin and the maxima

wx = 2 * pi * wx/H                             % Calculate wx and wy
wy = 2 * pi * wy/K

%% Computing X, the matrix containing I and A
A = ones(H * H,1);B = zeros(H);C = zeros(H); % creating empty matrixs                
for x=1:H
    for y=1:K
        pow = 1i * wx * x + 1i * wy * y;    % finding e^(i*wx*x+i*wy*y)
        B(x,y) = exp(pow);C(x,y) = exp(-pow);% ie matrix B and C
    end
end
I = double(reshape(moireImage,[],1));       % put the image into a vector
I=I./max(I);                                % normalize the image
M = [A,B(:),C(:)];                          % Put the vector in one matrix
X = M\I(:)                                 % Calculate Matrix X


%% Regenerate the original Moire pattern with matrix X
I2=X(1,1);  A2=2*X(2,1);                             % Get the parameters from matrix X
x=1:H; x=repmat(x',[1,K]); 
y=1:K; y=repmat(y,[H,1]);
moireImage2=I2+real(A2*exp(1i*wx*x+1i*wy*y));     % create the image
ft2 = fftshift(fft2(moireImage2));                  % Fourier Transform of Moire Pattern

%% Plotting the regenerated image and FT of it
    
figure();
subplot(2, 2, 1); imshow(abs(moireImage2));           % Moire Pattern

subplot(2,2,2); imshow(abs(ft2))                  % fft2 with fftshift

subplot(2,2,3); imshow(abs(ft2), [0, 300])        % Choosing a suitable value to show the maxima

subplot(2,2,4); imshow(abs(ft2), [0 10000])
