% Program to make Moire patterns with cos lines.
clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;

%% Encode the original image
I=imread("lena512.jpg");[H,K]=size(I);   % load and store image in the variable I
imshow(I)
I = double(I);
[y,x]=meshgrid(1:K,1:H);                                        % grid coordinates
w=2*pi/5;f=@(x)(1+tanh(5*sin(w*x)))/2;                          % periodic almost-rectangular function
theta=pi/3;moireImage=f(x*cos(theta)+y*sin(theta)+acos(I/255-0.5));     % phase-shifted grid image
subplot(3,3,1)
imshow(moireImage)

%% Computing X
[Imin,param]=imfrest2(moireImage,3,struct('numiter',10));        % Get the parameters of the image
for m=1:2
    for n=1:3
        if param.w(m,n) >= pi                                   % Get the correct values of wx and wy
            param.w(m,n) = param.w(m,n)-2*pi;
        end
    end
end 
X=[exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:))]\moireImage(:); % Compute X

%% Regenerate the original Moire pattern with matrix X
for n=1:3   
    if abs(imag(X(n,1))) < 10^(-5)
        I2 = real(X(n,1));                      % Decide which value is I and which is A
    end
    if imag(X(n,1)) > 10^(-10)
        A2 = 2*(X(2,1));
    end
end

% x=1:H; x=repmat(x',[1,K]); 
% y=1:K; y=repmat(y,[H,1]);
% moireImage2=I2+real(A2*exp(1i*(param.w(1,3))*x+1i*param.w(2,3)*y));     % The pure Moire patterns

f = @param.model;
moireImage2=f(x,y,param);     % The pure Moire patterns

%% Plotting the regenerated image and FT of it
    
subplot(3, 3, 2); imshow(moireImage2,[]); 
revealed_image = moireImage.*moireImage2;
subplot(3,3,3); imshow(revealed_image)          % Revealing original image

boxKernel = ones(10,10);                                    % Apply Simple Averaging filter
remove_moire = conv2(revealed_image, boxKernel, 'valid');
subplot(3,3,4);imshow(remove_moire,[])
title('Conv2 with 10x10 kernel')
%figure();imshow(remove_moire,[])

gauss = imgaussfilt(revealed_image,2.5,Padding="symmetric",FilterDomain='frequency');   % Try gaussian filter
remove_moire2 = imadjust(gauss);
subplot(3,3,5);imshow(remove_moire2,[])
title('imgaussfilt with sigma 2.5')

% %%
% ft = fftshift(fft2(revealed_image));                      % Try fft and ifft
% %figure()
% %imshow(ft)
% ft(abs(ft)>23) = 0;
% %remove_moire3 = ft.*lpf;
% %ft_moire = fftshift(fft2(moireImage2));
% remove_moire3 = ifft2(ifftshift(ft),'symmetric');
% %remove_moire3 = imadjust(remove_moire3);
% remove_moire3 = imadjust(remove_moire3);
% %subplot(2,3,6);imshow(remove_moire3)

N = 10;
remove_moire3 = imfilter(revealed_image,ones(N,N)/N^2);
remove_moire3 = imadjust(remove_moire3);
subplot(3,3,6);imshow(remove_moire3)
title('imfilter with 10x10 kernel')

filter = fspecial('gaussian',10,2.5);
remove_moire4 = imfilter(revealed_image,filter);
remove_moire4 = imadjust(remove_moire4);
subplot(3,3,7);imshow(remove_moire4,[])
title('imfilter with fspecial "gaussian" 10,2.5')
%%
