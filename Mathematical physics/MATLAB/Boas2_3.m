X=0;
Y=0;
[X,Y] = meshgrid(0:pi/100:pi);
Z = 0;
for k = 2:2:200
Z =Z + (k/(k^2-1))*exp(Y*(-k)).*sin(X*k);
end
Z = Z * (4/pi);
surf(X,Y,Z)