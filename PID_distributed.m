function [u,y] = PID_distributed(K, Ti, Td, sim_time, y_zad, u_prev, y_prev)

% Definicja sta³ych
global T;
global Umax ; 
global Umin;
Umax = 1 ; 
Umin = -1 ; 

%alokacja wektorów
y = ones(1,sim_time)*y_prev;
yzad = ones(1,sim_time)*y_zad;
ulocal = zeros(length(K),1); 
u = ones(1,sim_time)*u_prev;
sumU = 0 ; 
error = zeros(1,sim_time);
w = zeros(1,length(K));
name = "gbellmf"
%przeliczanie dla Pida
r2 = K.*Td./T;
r1 = K.*(T./(2.*Ti)-2.*Td./T-1);
r0 = K.*(1+T./(2.*Ti)+Td./T);

for i=3:sim_time
    %symulacja
    if i > 6 % bo wczeœniej nie mo¿na wyliczyæ wyjœcia, a trzeba liczyæ sterowanie
        y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
    end
    error(i) = yzad(i) - y(i)
    w = MembershipFun(name,Umin,Umax,length(K),u(i-1)); 
    for r = 1 : length(K)
        ulocal(r) = r2(r)*error(i-2)+r1(r)*error(i-1)+r0(r)*error(i) + u(i-1);
        sumU = sumU + ulocal(r)*w(r); 
    end
    sumU = sumU/sum(w); 
    u(i) = sumU;
    sumU = 0 ;
    %ograniczenie wartoœci
    if u(i) < Umin
        u(i) = Umin;
    elseif u(i) > Umax
        u(i) = Umax;
    end
   
end

end