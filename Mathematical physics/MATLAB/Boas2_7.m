[X,Y] = meshgrid(0:pi/100:pi,0:0.01:1);
Z = 0;
for k = 2:2:200
Z =Z + (k/((k^2-1)*sinh(k)))*sinh(k*(1-Y)).*sin(X*k);
end
Z = Z * (4/pi);
surf(X,Y,Z)