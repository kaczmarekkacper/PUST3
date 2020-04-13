function [u, y] = DMC(s, D, N, Nu, lambda, yzad, sim_time, u_prev, y_prev)

% Definicja sta³ych
global Umin;
global Umax;

%alokacja wektorów
y = ones(1,sim_time)*y_prev;
u = ones(1,sim_time)*u_prev;

%Obiekt DMC
dmc = DMCReg(s, D, N, Nu, lambda, Umin, Umax);
dmc.reset(u_prev);
dmc.setValue(yzad);

%Symulacja
for i=2:sim_time
    if i > 6
        y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
    end
    u(i) = dmc.countValue(y(i));
end

end

