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
%subplot(2,2,1)
imshow(moireImage)

%% Computing X
[Imin,param]=imfrest2(moireImage,3,struct('numiter',3));        % Get the parameters of the image
for m=1:2
    for n=1:3
        if param.w(m,n) >= pi                                   % Get the correct values of wx and wy
            param.w(m,n) = param.w(m,n)-2*pi;
        end
    end
end
%param.w
X=[exp(i*x(:)*param.w(1,:)+i*y(:)*param.w(2,:))]\moireImage(:); % Compute X

%% Regenerate the original Moire pattern with matrix X
for n=1:3   
    if abs(imag(X(n,1))) < 10^(-10)
        I2 = real(X(n,1));                      % Decide which value is I and which is A
    end
    if imag(X(n,1)) > 10^(-10)
        A2 = 2*(X(2,1));
    end
end

x=1:H; x=repmat(x',[1,K]); 
y=1:K; y=repmat(y,[H,1]);
moireImage2=I2+real(A2*exp(1i*(param.w(1,3))*x+1i*param.w(2,3)*y));     % The pure Moire patterns

%% Plotting the regenerated image and FT of it
    
% subplot(2, 2, 2); imshow(moireImage2);           
% subplot(2,1,2); imshow(moireImage.*moireImage2)          % Revealing original image

%% Animate the overlap process

for a=H:-5:1 
    %for b=1:K
        %imshow(moireImage.*padarray(moireImage2(a:H,1:K),[(a-1)/2]))
        
        imshow(moireImage.*padarray(moireImage2(a:H,1:K),[(a-1)],'post'))
    %end
end
for b=1:5:H
    imshow(moireImage.*padarray(moireImage2(b:H,1:K),[(b-1)],'post'))
end