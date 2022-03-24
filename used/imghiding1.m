%%
clc;
close all;  
clear;  
workspace;  
rows = 480;
columns = 480;
x0 = rows / 2;
y0 = columns / 2;
M = zeros(rows,columns);
n = 4;

%%
for x=1:rows
    for y=1:columns
%         if (x > x0||y <= y0)
%             angle = atan((y0 - y)/(x - x0));
%         end
%         if (x < x0||y <= y0)
%             angle = pi - atan((y0 - y)/(x - x0));
%         end
%         if (x < x0||y > y0)
%             angle = pi + atan((y0 - y)/(x - x0));
%         end
%         if (x > x0||y > y0)
%             angle = 2 * pi - atan((y0 - y)/(x - x0));
%         end
%         if (x == x0||y < y0)
%             angle = pi / 2;
%         end
%         if (x == x0||y > y0)
%             angle = 3 * pi / 2;
%         end
%         if (x == x0||y == y0)
%             angle = 0;
%         end
        if x < rows / 2
            lambda = 0.1;
            a = 0.0383;
        end
        if x >= rows / 2
            lambda = 0.18;
            a = 0.0689;  
        end
        k = 1;
        M(x,y) = 4 * (1 + cos(2 * pi / lambda * (atan((y0 - y)/(x - x0)) + a * sin(2 * pi * k / n)))) / 2 / n;
    end
end

imshow(M)

