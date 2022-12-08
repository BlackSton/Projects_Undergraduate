epsilon = 10;
sigma = 0.2;
r_c = sigma*2.5;
dr = r_c/1000;
r =0.1:dr:r_c;
V = 4*epsilon*((sigma./r).^12-(sigma./r).^6)-4*epsilon*((sigma/r_c)^12-(sigma/r_c)^6);
V_im = (24*epsilon/sigma)*((2*(-sigma./r).^13-(-sigma./r).^7)-(2*(-sigma/r_c)^13-(-sigma/r_c)^7));
plot(r,V_im)
hold on
plot(r,V)
xlim([0,r_c])
ylim([-30,25])
grid