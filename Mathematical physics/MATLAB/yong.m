[X,Y] = meshgrid(1:1:20);
Z = 0;
x1 = 4;
y1 = 10;
x2 = 10;
y2 = 10;
Z =Z + 1./sqrt((X-x1).^2+(Y-y1).^2);
Z =Z + 1./sqrt((X-x2).^2+(Y-y2).^2);
surf(X,Y,Z)