% Program to make Moire patterns with cos lines.
clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;
rows = 512;
columns = 512;

%% Creating a random Moire Pattern
I0=randn;  A=randn+1i*randn                     % generate random parameters
%I0=1;A=1+i
wx0=rand*2*pi-pi, wy0=rand*pi
H=512;K=512;                                    % size of image
x=1:H; x=repmat(x',[1,K]); 
y=1:K; y=repmat(y,[H,1]);
moireImage=I0+real(A*exp(1i*wx0*x+1i*wy0*y));   % create the image
ft = fftshift(fft2(moireImage));                % Fourier Transform of Moire Pattern

%% Plotting the original image and FT of it

subplot(2, 2, 1);imshow(moireImage);

subplot(2,2,2);imshow(abs(ft))

subplot(2,2,3);imshow(abs(ft), [0, 300])

subplot(2,2,4);imshow(abs(ft), [0 1000])

%% Finding the location of the maximum point
% Set the magnitude of origin to 0
% J = abs(ft)/1000;
% origin = rows/2+1;
% J(origin,origin) = 0;
% 
% [jx, jy] = find(ismember(J, max(J(:))));        % locate the maxima of the ft.
% 
% wx = jx(2) - origin;                            % find the difference between
% wy = jy(2) - origin;                            % the origin and the maxima
% 
% wx = 2 * pi * wx/H                             % Calculate wx and wy
% wy = 2 * pi * wy/K

%% Computing X
% A = ones(H * K,1);                                 % creating empty matrixs
% B = zeros(H);C = zeros(H);
% for xx=1:H
%     for yy=1:K
%         pow = 1i * wx * xx + 1i * wy * yy;           % finding e^(i*wx*x+i*wy*y)
%         B(xx,yy) = exp(pow);                          % ie matrix B and C
%         C(xx,yy) = exp(-pow);
%     end
% end
% 
% % Put the matrix in a vector
% I = reshape(moireImage,H*K,1);                      % put the image into a vector
% M = [A,B(:),C(:)];                                  % Put the vector in one matrix       
% X = M\I

%X=[ones(K*H,1),exp(i*wx0*x(:)+i*wy0*y(:)),exp(-i*wx0*x(:)-i*wy0*y(:))]\moireImage(:)
%=[ones(K*H,1),exp(i*wx*x(:)+i*wy*y(:)),exp(-i*wx*x(:)-i*wy*y(:))]\moireImage(:)
[Imin,param]=imfrest2(moireImage,3,struct('numiter',3));
X=[exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:))]\moireImage(:)

% %% Regenerate the original Moire pattern with matrix X
I2=real(X(1,1));  A2=2*X(2,1)                             % Get the parameters from matrix X
x=1:H; x=repmat(x',[1,K]); 
y=1:K; y=repmat(y,[H,1]);
moireImage2=I2+real(A2*exp(1i*(param.w(1,2)-2*pi)*x+1i*param.w(2,2)*y));     % create the image
ft2 = fftshift(fft2(moireImage2));                  % Fourier Transform of Moire Pattern

%% Plotting the regenerated image and FT of it
    
figure();
subplot(2, 2, 1); imshow(moireImage2);           % Moire Pattern

subplot(2,2,2); imshow(abs(ft2))                  % fft2 with fftshift

subplot(2,2,3); imshow(abs(ft2), [0, 300])        % Choosing a suitable value to show the maxima

subplot(2,2,4); imshow(abs(ft2), [0 1000])