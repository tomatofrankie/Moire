clear all

% Plotting vertical lines
x = 1:100;
figure()
plot([x; x], [max(x) * zeros(1, length(x))*min(ylim); max(x) * ones(1, length(x))*max(ylim)], "Color", [0,0,0], "LineWidth", 1)

% Plotting the inclined lines
slope = 1;
%slope = tan((90 - 45)/180 * pi)
min = -100 * abs(slope - 1) + 1;
for con=min:100
    k = slope * x + con;
    %y = slope * x + con;
    y(k > max(x)) = max(x);
    y(k < 0) = 0;
    hold on
    plot(x, y, "LineWidth", 1, "Color", [0,0,0])
end

%% 
clear all

% Plotting vertical lines
x = 1:100;
figure()
plot([x; x], [max(x) * zeros(1, length(x))*min(ylim); max(x) * ones(1, length(x))*max(ylim)], "Color", [0,0,0], "LineWidth", 1)

% Plotting the inclined lines
slope = 10;
angle = atand(slope);
%slope = tan(deg2rad(90 - angle))
%min = -100 * abs(slope - 1) + 1
c = 1 / cos(deg2rad(angle))
for con=-slope * 100:100
    k = slope * x + c * con;
    y = slope * x + c * con;
    y(k > max(x)) = max(x);
    y(k < 0) = 0;
    hold on
    plot(x, y, "LineWidth", 1, "Color", [0,0,0])
end
    