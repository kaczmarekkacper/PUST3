function [u, y] = DMC_distributed(s, D, N, Nu, lambda, yzad, sim_time, u_prev, y_prev)

% Definicja sta�ych
global Umin;
global Umax;
global Rmax; 
global Rmin; 
global memFun; 
global Rvar; 

%alokacja wektor�w
y = ones(1,sim_time)*y_prev;
u = ones(1,sim_time)*u_prev;
ulocal = zeros(1,length(D));
sumU = 0 ; 
%Obiekt DMC
for n = 1 : length(D) 
    dmcloc(n) = DMCReg(s(1:end,n), D(n), N(n), Nu(n), lambda, -10e19, 10e19);
    dmcloc(n).reset(u_prev);
    dmcloc(n).setValue(yzad);
end

%Symulacja
for i=2:sim_time
    if i > 6
        y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
    end
    if strcmp(Rvar,'u(i-1)')
        w = MembershipFun(memFun,Rmin,Rmax,length(D),u(i-1));
    else
        w = MembershipFun(memFun,Rmin,Rmax,length(D),y(i));
    end
    for n = 1 : length(D)
        ulocal(1,n) = dmcloc(n).countValue(y(i));
        sumU = ulocal(1,n)*w(n) + sumU ;
    end
    
    u(i) = sumU/sum(w);
%     sum(w)
%     sumU
%     u(i)
    sumU = 0 ;
   %ograniczenie warto�ci
    if u(i) < Umin
        u(i) = Umin;
    elseif u(i) > Umax
        u(i) = Umax;
    end
end

end

