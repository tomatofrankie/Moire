M = 512; N = M;
I = imread('lena512.jpg');
I = imresize(I, [M N]);
x = 1:M;
x = x' * ones(1, N);
y = 1:N;
y = ones(M, 1) * y;
w = 2 * pi/5;
f = @(x)tanh(sin(w * x));
% f = @(x)sin(w * x)
theta = pi/3;
I1 = f(x * cos(theta) + y * sin(theta) + acos(double(I)/255 - 0.5));
theta = pi/3;
I2 = f(x * cos(theta) + y * sin(theta));
s = 1;
%exp(-s^2 * (wx .^ 2 + wy .^ 2)/2);
I1_1 = I1 .* I2;
imshow(I1_1)
imcontrast()
filtered = imfilter(I1 .* I2, exp(-s^2 * (x .^ 2 + y .^ 2)/2));

%imcontrast(I1, I2))