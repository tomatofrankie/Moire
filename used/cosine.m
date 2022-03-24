% Program to make Moire patterns with cos lines.

clc;    
close all;  
clear;  
workspace;  
format long g;
format compact;
fontSize = 24;
rows = 480;
columns = 640;

% Make an image with rippled lines in it.
rowVector = (1 : rows)';
period = 20; 
amplitude = 0.6; 
offset = 1 - amplitude; 
cosVector = amplitude * (1 + cos(2 * pi * rowVector / period))/2 + offset;
ripplesImage = repmat(cosVector, [1, columns]);
subplot(1, 3, 1);
imshow(ripplesImage, []);
axis on;
title('Ripple image', 'FontSize', fontSize);

% Make an image with tilted rippled lines in it.
ripplesImage2 = imrotate(ripplesImage, 10, 'crop'); 
subplot(1, 3, 2);
imshow(ripplesImage2, []);
axis on;
title('Tilted Ripples', 'FontSize', fontSize);

% Multiply the ripple images together to get a Moire pattern.
grayImage = ripplesImage .* double(ripplesImage2);
subplot(1, 3, 3);
imshow(grayImage, []);
axis on;
title('Product = Moire Image', 'FontSize', fontSize);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); 
drawnow;