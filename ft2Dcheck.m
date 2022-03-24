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
I0=rand,  A=rand+1i*rand                    % generate random parameters
% I0=1;A=1+i
wx0=rand*2*pi-pi, wy0=2*rand*pi-pi
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


%% Computing X
[Imin,param]=imfrest2(moireImage,3,struct('numiter',3));
for m=1:2
    for n=1:3
        if param.w(m,n) >= pi
            param.w(m,n) = param.w(m,n)-2*pi;
        end
    end
end
param.w
X=[exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:))]\moireImage(:)

% %% Regenerate the original Moire pattern with matrix X
for n=1:3
    if abs(imag(X(n,1))) < 10^(-10)
        I2 = real(X(n,1));
    end
    if imag(X(n,1)) > 10^(-10)
        A2 = 2*(X(2,1));
    end
end

%A2=2*X(2,1)                             % Get the parameters from matrix X
x=1:H; x=repmat(x',[1,K]); 
y=1:K; y=repmat(y,[H,1]);
moireImage2=I2+real(A2*exp(1i*(param.w(1,2))*x+1i*param.w(2,2)*y));     % create the image
ft2 = fftshift(fft2(moireImage2));                  % Fourier Transform of Moire Pattern

%% Plotting the regenerated image and FT of it
    
figure();
subplot(2, 2, 1); imshow(moireImage2);           % Moire Pattern

subplot(2,2,2); imshow(abs(ft2))                  % fft2 with fftshift

subplot(2,2,3); imshow(abs(ft2), [0, 300])        % Choosing a suitable value to show the maxima

subplot(2,2,4); imshow(abs(ft2), [0 1000])