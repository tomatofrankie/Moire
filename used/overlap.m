clc
%close all
clear 

figure(Name='')
original = imread('lena512.jpg');
subplot(1,3,1)
imshow(original, [])

rotate = imrotate(original,5,'crop');
subplot(1,3,2)
imshow(rotate, [])

overlaped1 = original .* rotate;
subplot(1,3,3)
imshow(overlaped1, [])
