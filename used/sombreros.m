% Program to make Moire patterns.
% First Moire pattern is from interference (multipication) to two linear sine wave ripples.
% This is equivalent to the Fourier diffraction pattern of two infinitely long slits.
% The second Moire pattern is from interference (multipication) to two Sombrero function ripples.
% This is equivalent to the Fourier diffraction pattern of two infinitely small circular apertures.

clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format long g;
format compact;
fontSize = 24;
rows = 480;
columns = 640;

% Make an image with rippled lines in it.
rowVector = (1 : rows)';
period = 20; % 20 rows
amplitude = 0.5; % Magnitude of the ripples.
offset = 1 - amplitude; % How much the cosine is raised above 0.
cosVector = amplitude * (1 + cos(2 * pi * rowVector / period))/2 + offset;
ripplesImage = repmat(cosVector, [1, columns]);
subplot(2, 3, 1);
imshow(ripplesImage, []);
axis on;
title('Ripple image', 'FontSize', fontSize);

% Make an image with tilted rippled lines in it.
ripplesImage2 = imrotate(ripplesImage, 10, 'crop'); % Rotate first image.
subplot(2, 3, 2);
imshow(ripplesImage2, []);
axis on;
title('Tilted Ripples', 'FontSize', fontSize);

% Multiply the ripple images together to get a Moire pattern.
grayImage = ripplesImage .* double(ripplesImage2);
subplot(2, 3, 3);
imshow(grayImage, []);
axis on;
title('Product = Moire Image', 'FontSize', fontSize);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]); % Maximize figure.
drawnow;

% Get a sombrero function.  http://en.wikipedia.org/wiki/Sombrero_function
% This is the diffraction pattern of a circular aperture.
xo = (1 : columns) - columns/2;
yo = (1 : rows) - rows/2;
% Scale pattern
scaleFactor = 50;
x = xo / scaleFactor;
y = yo / scaleFactor;
sombrero = zeros(rows, columns);
[X, Y] = meshgrid(x, y);
radii = sqrt(X.^2 + Y.^2);
sombrero = 2 * besselj(1, pi*radii) ./ (pi * radii);
subplot(2, 3, 4);
imshow(sombrero, []);
axis on;
title('Sombrero Function', 'FontSize', fontSize);

% Now, same thing but shift the center of the sombrero
xo = (1 : columns) - 1.2 * columns/2;
yo = (1 : rows) - 1.5 * rows/2;
% Scale pattern
x = xo / scaleFactor;
y = yo / scaleFactor;
shiftedSombrero = zeros(rows, columns);
[X, Y] = meshgrid(x, y);
radii = sqrt(X.^2 + Y.^2);
shiftedSombrero = 2 * besselj(1, pi*radii) ./ (pi * radii);
subplot(2, 3, 5);
imshow(shiftedSombrero, []);
axis on;
title('Shifted Sombrero Function', 'FontSize', fontSize);

% Multiply the ripple images together to get a Moire pattern.
grayImage2 = sombrero .* double(shiftedSombrero);
subplot(2, 3, 6);
imshow(grayImage2, []);
axis on;
title('Product = Moire Image', 'FontSize', fontSize);

