[X,Y] = meshgrid(0:0.1:10);
Z = 0;
for k = 2:2:200
Z =Z + (1/k)*exp(Y*-(k*pi*1/10)^2).*sin(X*(k*pi/10));
end
Z = 100-10*X-Z*(400/pi);
surf(X,Y,Z)