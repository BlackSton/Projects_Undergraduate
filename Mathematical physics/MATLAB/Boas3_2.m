v = VideoWriter('Boas3_2.avi');
open(v);
x = 0:0.01:10; % (0,10) 영역안에서
for t = 0:0.1:20 % 10초까지
    y_t = 0;
    for k = 1:2:200
        y_t =y_t + (1/k)*exp(t*-(k*pi*1/10)^2).*sin(x*(k*pi/10));
    end
    y = y_t*(400/pi);
 plot(x,y);
 ylim([0 100])
 F = getframe(gca);
writeVideo(v,F);
end
 close(v);