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
moireImage0=f(x*cos(theta)+y*sin(theta));
subplot(3,3,1)
imshow(moireImage)

%% Regenerate the original Moire pattern with matrix X

[Imin,param]=imfrest2(moireImage,3,struct('numiter',10));
for m=1:2
    for n=1:3
        if param.w(m,n) >= pi                                   % Get the correct values of wx and wy
            param.w(m,n) = param.w(m,n)-2*pi;
        end
    end
end 
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
filter = fspecial('laplacian',1);
remove_moire4 = imfilter(revealed_image,filter);
% remove_moire4 = medfilt2(revealed_image,[8,8]);
remove_moire4 = imadjust(remove_moire4);
imshow(remove_moire4,[])

%% Try removing spikes in fft
ft = fftshift(fft2(revealed_image));                      % Try fft and ifft
% figure()
% subplot(2,2,1),imshow(abs(ft),[0,10000])

jx = param.w(1,3)*H/(2*pi);jy = param.w(2,3)*K/(2*pi);      % Calculate location of maximum point
centrex = H/2;centrey = K/2;
% jx1 = round(centrex + jx);jx2 = round(centrex - jx);
% jy1 = round(centrey + jy);jy2 = round(centrey - jy);

% Set maximum points to zero

ft(1:centrex+round(jx/2),1:centrey+round(jy/2)) = 0;ft(centrex-round(jx/2):H,centrey-round(jy/2):K) = 0;

% ft(:,1:centrey+round(jy/2)) = 0;ft(:,centrey-round(jy/2):K) = 0;
% ft(1:centrex+round(jx/2),:) = 0;ft(centrex-round(jx/2):H,:)=0;

% subplot(2,2,2),imshow(abs(ft),[0,10000])
remove_moire5 = abs(ifft2(ifftshift(ft)));
remove_moire5 = imadjust(remove_moire5);
subplot(338);imshow(remove_moire5)
title('removing spikes from ft')
subplot(339),imshow(imadjust(imgaussfilt(remove_moire5,1.5,Padding="symmetric",FilterDomain='frequency')),[])
title('imgaussfilt after removing spikes')