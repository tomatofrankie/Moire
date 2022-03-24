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
moireImage2 = cell(ceil((D+1)/H),ceil((E+1)/K));
gridx = 1;
[y,x]=meshgrid(1:K,1:H);
moireImage=double(rgb2gray(moireImage0));                   % convert to grayscale image
for a=1:H:D                                      % move the kenerl
    gridy = 1;
    for b=1:K:E
        moireImage=imcrop(moireImage,[b a H-1 K-1]);        % crop the image to the size of kernel
        if size(moireImage,1) < H-1
            moireImage = padarray(moireImage,H-size(moireImage,1),1,'post');
        end
        if size(moireImage,2) < K-1
            moireImage = padarray(moireImage,[0,K-size(moireImage,2)],1,'post');
        end

        %% Finding the location of the maximum point
        [Imin,param]=imfrest2(moireImage,3,struct('numiter',10));        % Get the parameters of the image
        for m=1:2
            for n=1:3
                if param.w(m,n) >= pi                                   % Get the correct values of wx and wy
                    param.w(m,n) = param.w(m,n)-2*pi;
                end
            end
        end
        %% Regenerate the original Moire pattern with matrix X
        for n=1:3   
            if abs(imag(X(n,1))) < 10^(-1)
                I2 = real(X(n,1));                      % Decide which value is I and which is A
            end
            if imag(X(n,1)) > 10^(-1)
                A2 = 2*(X(2,1));
            end
        end
        
        x=1:H; x=repmat(x',[1,K]); 
        y=1:K; y=repmat(y,[H,1]);
        moireImage2{gridx,gridy}=I2+real(A2*exp(1i*(param.w(1,3))*x+1i*param.w(2,3)*y))     % The pure Moire patterns
        gridy = gridy + 1;
    end
    gridx = gridx + 1;
end