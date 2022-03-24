% Program to get the X matrix of an image
clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%% Load an image
moireImage0= imread("moire1.jpg");      % read an image            
[D,E,~]=size(moireImage0);              % get the size of the image

%% Crop the image into small squares
H=50;K=50;                                  % size of kernel
for a=1:D-H+1                               % move the kenerl
    for b=1:E-K+1
        moireImage=rgb2gray(moireImage0);                   % convert to grayscale image
        moireImage=imcrop(moireImage,[b a H-1 K-1]);        % crop the image to the size of kernel
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
        J(xorigin+1,yorigin) = 0;           % Set the magnitude of origin to 0
        
        [J0, jx] = max(J);                  % locate the maxima of the ft.
        [J1, jy] = max(pagetranspose(J));        
        [~,jx] = max(J1);[~,jy] = max(J0);
       
        wx = xorigin - jy;wy = yorigin - jx;   % find the difference between the origin and the maxima
        wx = 2 * pi * wx/H;wy = 2 * pi * wy/H; % calculate wx and wy
        
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
        X = M\I(:);                                 % Calculate Matrix X
    end
end