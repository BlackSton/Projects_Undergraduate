x=1:0.1:19.9;
y=1:1:18;
[X,Y] = meshgrid(x,y);
Z = csvread('output.csv',1,0);
surf(X,Y,Z)