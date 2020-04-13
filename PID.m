function [u,y] = PID(K, Ti, Td, sim_time, y_zad, u_prev, y_prev)

% Definicja sta³ych
global Umin;
global Umax;
global T;
global Umaxchange;

%alokacja wektorów
y = ones(1,sim_time)*y_prev;
yzad = ones(1,sim_time)*y_zad;
u = ones(1,sim_time)*u_prev;
error = zeros(1,sim_time);

%przeliczanie dla Pida
r2 = K*Td/T;
r1 = K*(T/(2*Ti)-2*Td/T-1);
r0 = K*(1+T/(2*Ti)+Td/T);

for i=3:sim_time
    %symulacja
    if i > 6 % bo wczeœniej nie mo¿na wyliczyæ wyjœcia, a trzeba liczyæ sterowanie
        y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
    end
    error(i) = yzad(i) - y(i);
    u(i) = r2*error(i-2)+r1*error(i-1)+r0*error(i) + u(i-1);
    
    du = u(i)-u(i-1);
    
    %ograniczenie zmiany sterownaia
    if du > Umaxchange
        u(i) = u(i-1) + Umaxchange;
    elseif du < -Umaxchange
        u(i) = u(i-1) - Umaxchange;
    end
    %ograniczenie wartoœci
    if u(i) < Umin
        u(i) = Umin;
    elseif u(i) > Umax
        u(i) = Umax;
    end
end

end

