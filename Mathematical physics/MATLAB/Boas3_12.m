v = VideoWriter('Boas3_12.avi');
open(v);
x = 0:0.01:1; % (0,1) �����ȿ���
for t = 0:0.01:20 % 10�ʱ���
    y_t = 0;
    for k = 1:2:1000
        y_t= y_t+(1/(k*(4-k^2)))*sin(k*pi*x).*exp(-1j*(k^2)*t);
    end
    y = y_t*(8/pi);
    y = abs(y).^2;
 plot(x,y);
 F = getframe(gca);
writeVideo(v,F);
end
 close(v);