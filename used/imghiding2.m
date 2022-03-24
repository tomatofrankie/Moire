clc;
close all;  
clear;  
workspace;  
rows = 480;
columns = 480;
x0 = rows / 2;
y0 = columns / 2;
M = zeros(rows,columns);
F = zeros(rows,columns);

for x=1:rows
    for y=1:columns
        if x < rows / 2
            lambda1 = 0.1;
            a1 = 0.0383;
        end
        if x >= rows / 2
            lambda1 = 0.18;
            a1 = 0.0689;  
        end
        k = 0;
        M(x,y) = (1 + cos(2 * pi / lambda1 * (atan((y0 - y)/(x - x0)) ...
            + a1 * sin(2 * pi * k / 1)))) / 2;
    end
end
imshow(M)
