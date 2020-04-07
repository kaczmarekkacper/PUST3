function [u, y] = DMC(s, D, N, Nu, lambda, yzad, sim_time, u_prev, y_prev)

% Definicja sta³ych
global Upp;
global Ypp;
global Umin;
global Umax;
global Umaxchange;
Umaxchange = 0.2;

%alokacja wektorów
y = ones(1,sim_time)*y_prev;
u = ones(1,sim_time)*u_prev;

%Obiekt DMC
dmc = DMCReg(s, D, N, Nu, lambda, Umin, Umax, Umaxchange);
dmc.reset(u_prev);
dmc.setValue(yzad);

%Symulacja
for i=2:sim_time
    if i > 12
        y(i) = symulacja_obiektu7Y(u(i-10),u(i-11),y(i-1),y(i-2));
    end
    u(i) = dmc.countValue(y(i));
end

end

