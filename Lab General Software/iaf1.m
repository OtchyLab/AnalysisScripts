function iaf1(Ie,Vth,dt,T)

EL = -65;
Vreset = -65;
Vth = -50;
taum = 10;
Rm = 10;	%M Ohm
s = taum/dt;
v = EL;
vshow(1) = v;
t = 0;
j = 1;

while t < T

    j = j + 1;
    t = (j-1)*dt;
    I = 0;
    if t > 1,
       I = Ie;		% nanoAmps
    end
    if v < Vth
       v = (s*v+EL+Rm*I)/(s+1);
       vshow(j) = v;
    else
       v = Vreset;
       vshow(j) = 10;
    end

end

tim = 0:dt:(j-1)*dt;
plot(tim,vshow)
xlabel('t (ms)','fontsize',16)
ylabel('V (mV)','fontsize',16)
head = ['iaf1(' num2str(Ie) ',' num2str(Vth) ','];
head = [head num2str(dt) ',' num2str(T) ')'];
title(head,'fontsize',16)