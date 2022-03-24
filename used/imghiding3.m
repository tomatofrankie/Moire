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
F = zeros(rows,columns);
n = 50;

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
            lambda1 = 0.1;
            a1 = 0.0383;
        end
        if x >= rows / 2
            lambda1 = 0.18;
            a1 = 0.0689;  
        end

        k = 0;

        M(x,y) = (1 + cos(2 * pi / lambda1 * (atan((y0 - y)/(x - x0)) + a1 * sin(2 * pi * k / 1)))) / 2;
        %F(x,y) = (1 + cos(2 * pi / lambda2 * (atan((y0 - y)/(x - x0)) + a2 * sin(2 * pi * k / n)))) / 2;
        
        if (x > x0||y <= y0)|(x < x0||y <= y0)|(x < x0||y > y0)|(x > x0||y > y0)
            angle = atan((y0 - y)/(x - x0));
        elseif (x == x0||y < y0)
            angle = pi / 2;
        elseif (x == x0||y > y0)
            angle = 3 * pi / 2;
        elseif (x == x0||y == y0)
            angle = 0;
        else
            
        end
%%
        if x < rows / 2
            F(x,y) = M(x,y);
            for k=0:n-1
                %F(x,y) = M(x,y);
                F(x,y) = F(x,y) + (1 + cos(2 * pi / lambda1 * (angle + a1 * sin(2 * pi * k / n)))) / 2;
            end
            F(x,y) = F(x,y) / n;
        else
            F(x,y) = M(x,y);
        end
%%
        if x >= rows / 2
            C(x,y) = M(x,y);
            for k=0:n-1
                %F(x,y) = M(x,y);
                C(x,y) = C(x,y) + (1 + cos(2 * pi / lambda1 * (angle + a1 * sin(2 * pi * k / n)))) / 2;
            end
            C(x,y) = C(x,y) / n;
        else
            C(x,y) = M(x,y);
        end

%%
        D(x,y) = M(x,y);
        %D(x,y) = 0;
        for k=0:n-1
            %D(x,y) = M(x,y);
            D(x,y) = D(x,y) + M(x,y);
        end
        D(x,y) = D(x,y) / n / 2;
        
    end
end

%%
subplot(2,2,1)
imshow(M)
%%
subplot(2,2,2)
imshow(F)
%%
subplot(2,2,3)
imshow(C)
%%
subplot(2,2,4)
imshow(D)

%imshow(F)