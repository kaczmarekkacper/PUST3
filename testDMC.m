function [] = testDMC(nr,lambda)
figures = 1 ; 
global Ypp ;
global k ; 
global Upp; 
global Ymin ; 
global Ymax ;
global funcall ;
funcall = 1 ;
tau = 10 ; 
y = ones(k, nr)*Ypp; %alokacja wektora o d�ugo�ci symulacji
u = ones(k,nr)*Upp; %Sterowanie pocz�tkowe 
yDMC = zeros(k,nr); 
uDMC = zeros(k,nr); 
s = zeros(k,nr); 
rng(19981998); % ustalenie seedu, żeby rand zawsze dawał to samo przy każdym uruchomieniu

for i = 0:9
    y_zad_traj(1+200*i:200*(i+1)) = Ymin + (Ymax-Ymin)*rand(1);
end
%ustawianie "centr�w" funkcji przynalezno�ci regulator�w przy rozk�adzie
%r�wnomiernym 
%centerDistance = (Rmax-Rmin)/(nr-1); 
%centers = Rmin : centerDistance : Rmax; 

%r�czne ustawianie centr�w funkcji przynale�no�ci
if(nr == 2)
   centers = [-1;1];
elseif( nr == 3)
   centers = [-1;-0.4;0.9]; 
elseif(nr == 4)
   centers = [-1;-0.4;0.3;0.9];
elseif(nr == 5)
   centers = [-1;-0.55;0;0.2;0.9];
end

%Zbieranie lokalnych odpowiedzi skokowych obiektu
ypp = 0 ; %warto�� pocz�tkowa wyj�cia przy danym skoku
upp = 0 ; %warto�c pocz�tkowa sterowania przy danym skoku
skok = 0.05; %warto�� skoku

for reg = 1 : nr
    %wyznaczenie punkt�w pracy w centrach funkcji przynale�no�ci
    u(10:k,reg) = centers(reg);
    for i = 7:k
        y(i,reg) = symulacja_obiektu7y(u(i-5,reg),u(i-6,reg),y(i-1,reg),y(i-2,reg));
    end
    upp = centers(reg) ;
    ypp = y(i,reg);
    
    %ustalenie warto�ci pocz�tkowych w punkcie pracy
    u(1:10,reg) = ones(10,1)*upp;
    y(1:end,reg) = ones(k,1)*ypp;
    if upp == 1 
        skok = -1*skok; 
    end
    %ma�y skok w celu zebrania lokalnej odpowiedzi    
    u(10:k,reg) = upp + skok; 
    for i = 7:k
        y(i,reg) = symulacja_obiektu7y(u(i-5,reg),u(i-6,reg),y(i-1,reg),y(i-2,reg));
    end
    
    %przeskalowanie aby otrzyma� odpowied� na skok jednostkowy
    yDMC(1:end-tau,reg) = y(tau+1:end,reg); 
    yDMC(end-tau:end,reg) = y(end-tau:end,reg);
    yDMC(1:end,reg) = yDMC(1:end,reg) - ypp;
    uDMC(1:end-tau,reg) = u(tau+1:end,reg);
    uDMC(end-tau:end,reg) = u(end-tau:end,reg);
    uDMC(1:end,reg) = uDMC(1:end,reg) - upp;
    s(1:end,reg) = yDMC(1:end,reg)/skok;
%     plot(yDMC(1:end,reg))
%     hold on
%     plot(uDMC(1:end,reg))
%     hold on
%     figure
%     plot(s(1:end,reg))
%     hold on
%     if centers(reg) ~= 0
%         plot(uDMC(1:end,reg)/skok)
%     end
end


%Parametry regulator�w 
D = setD(s);
N = D;
Nu = fix(N/2);

%Symulacja procesu
for i = 0:9
    yzad = y_zad_traj(1+200*i);
    if i == 0
        [u, y] = DMC_distributed(s, D, N, Nu, lambda, yzad, 200, Upp,Ypp);
    else
        [u, y] = DMC_distributed(s, D, N, Nu, lambda, yzad, 200, u_traj(200*i), y_traj(200*i));
    end
    u_traj(1+200*i:200*(i+1)) = u;
    y_traj(1+200*i:200*(i+1)) = y;
end

%B��d �redniokwadratowy
E = (y_traj-y_zad_traj)*(y_traj-y_zad_traj)';

%Rysowanie wynik�w symulacji
k = 2000;
title = sprintf('DMC N=D i Nu=N/2 lambda=[%.2f,%.2f] E=%.2f ',lambda(1),lambda(2),E);
if figures
    plotProcess(u_traj, y_traj,title);
end
subplot(2,1,1);
hold on;
plot(1:k, y_zad_traj, 'r:');
legend('Wyj�cie procesu', 'Warto�� zadana', 'Location', 'southwest');
hold off;
k = 250;
%Zapisanie wynik�w
% mkdir('results/5/DMC')
matlab2tikz(sprintf('results//5//DMC//DMC_reg%d_N=D_Nu=0,5N_lambda=%.2f_%.2f.tex',nr,lambda(1),lambda(2)))


end

