syms x y;
x1 = 0; %xmin
x2 = 1 - y/2; %xmax
y1=0;%ymin
y2=2;%ymax
f = 1 + y; % function
result1 = int(f,[x1 x2]); %first integral
result2 = int(result1,[y1 y2]); %second integral
result2