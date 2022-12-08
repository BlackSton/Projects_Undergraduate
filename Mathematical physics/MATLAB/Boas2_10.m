[X,Y] = meshgrid(0:0.1:10);
Z = 0;
for k = 1:2:200
Z =Z + (400/(k*pi*sinh(k*pi)))*(sinh((k*pi/10)*(10-Y)).*sin(X*(k*pi/10))+sinh((k*pi/10)*(10-X)).*sin(Y*(k*pi/10)));
end
Z(51,51) % (5,5)
surf(X,Y,Z)