v = VideoWriter('Boas3_11.avi');
open(v);
x = 0:pi/100:pi; % (0,pi) 영역안에서
for t = 0:0.01:10 % 10초까지
    y_t = 0;
    for k = 1:1000
        y_t= y_t+(1/k)*sin(k*x).*exp(-1j*(k^2)*t);
    end
    y = y_t*(4/pi);
    y = abs(y).^2;
 plot(x,y);
 F = getframe(gca);
writeVideo(v,F);
end
 close(v);