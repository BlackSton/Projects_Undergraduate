L = 1:1:20;
v = 1:1:20;
m = 1;
% 시간 설정 기준 Time = L*1.5/v

[X,Y] = meshgrid(L,v);
Z = 0;
for k = 2:2:200
Z =Z + (k/((k^2-1)*sinh(k)))*sinh(k*(1-Y)).*sin(X*k);
end
Z = Z * (4/pi);
surf(X,Y,Z)
