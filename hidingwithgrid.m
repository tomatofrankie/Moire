
I=imread("lena512.jpg");[M,N]=size(I);   % load and store image in the variable I
I = double(I);
[y,x]=meshgrid(1:N,1:M);                                        % grid coordinates
w=2*pi/5;f=@(x)(1+tanh(5*sin(w*x)))/2;                          % periodic almost-rectangular function
theta=pi/3;I1=f(x*cos(theta)+y*sin(theta)+acos(I/255-0.5));     % phase-shifted grid image
theta=pi/3;I2=f(x*cos(theta)+y*sin(theta));                     % pure grid image


%figure(1),image(255*I1),colormap(gray(256)),axis image,axis off,title('Phase-shifted grid image')
%figure(2),image(255*I2),colormap(gray(256)),axis image,axis off,title('Pure grid image')
%figure(3),image(255*I1.*I2),colormap(gray(256)),axis image,axis off,title('Moiré effect (multiplication of the grid images)')

figure()
subplot(2,2,1)
imshow(255*I1)
subplot(2,2,2)
imshow(255*I2)
subplot(2,2,3)
imshow(255*I1.*I2)