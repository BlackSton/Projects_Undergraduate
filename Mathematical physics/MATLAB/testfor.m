v = VideoWriter('peaks.avi');
open(v);
x = 0:pi/100:pi;
h = 1.054571800*10^(-34); %디랙 상수
m = 0.001/(6.02*10^23); % 수소 원자의 질량
for t = 0:0.1:10
    for k = 1:2:100
        y= (1/k)*sin(k*x).*exp(-1i*(k^2)*t);
    end
    y = y*(4/pi);
    y = abs(y).^2;
 plot(x,y);
 F = getframe(gca);
writeVideo(v,F);
end
 close(v);