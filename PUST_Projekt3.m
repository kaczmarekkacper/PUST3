%% 0. Definicje makr 

clear;
close all;
more off;

global Upp;
global Ypp;
global Umin;
global Umax;
global Ymax;
global Ymin;
global T;
global k; 
%zmienna decyzyjna, prze��czaj�ca regulatory
global Rvar;
%zakres warto�ci zmiennej decyzyjnej
global Rmax; 
global Rmin; 
%funkcja przynale�no�ci
global memFun; 



Upp = 0;
Ypp = 0;
Umin = -1;
Umax = 1;
Ymax = 0.149670193432045; % wyznaczone eksperymentalnie
Ymin = -4.774746119647459;
T = 0.5;
k = 250;

figures = 1; % czy wyświetlać wykresy
saving = 1; % czy zapisywać dane

%% 1. Sprawdzic poprawnosc wartosci Upp, Ypp.

mkdir('results/1');

y = zeros(k, 1); %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
for i = 3:k
    y(i) = symulacja_obiektu7y(Upp,Upp,y(i-1),y(i-2));
end

%rysowanie wykresów
if figures
    plotProcess(u, y,'Sprawdzenie punktu pracy');
end

%zapisanie plików do tex, aby wrzucić wykres w tikz do latex
if saving
    matlab2tikz('results/1/CheckWorkpoint.tex');
end

%% 2. Wyznaczyc symulacyjnie odpowiedzi skokowe procesu dla kilku zmian sygnału sterujacego,
%     przy uwzglednieniu ograniczen wartosci tego sygnału, jego wartosc na poczatku
%     eksperymentu wynosi 0. Narysowac te odpowiedzi na jednym rysunku. Narysowac
%     charakterystyke statyczna procesu y(u). Czy własciwosci statyczne i dynamiczne procesu
%     sa (w przyblizeniu) liniowe?

mkdir('results/2');

% 2.1 Odpowiedzi skokowe

%Skok z 0 na 1 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = 1;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    figure;
    subplot(2, 1, 1);
    stairs(1:k, y);
    title('Odpowiedzi skokowe');
    xlabel('k');
    ylabel('y');
    hold on;
    subplot(2, 1, 2);
    stairs(1:k, u);
    xlabel('k');
    ylabel('u');
    hold on;
end

%Skok z 0 na 0.5 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = 0.5;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    subplot(2, 1, 1);
    stairs(1:k, y);
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u);
    xlabel('k');
    ylabel('u');
end

%Skok z 0.0 na 0.2 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = 0.2;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    subplot(2, 1, 1);
    stairs(1:k, y);
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u);
    xlabel('k');
    ylabel('u');
end

%Skok z 0.0 na -0.2 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = -0.2;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    subplot(2, 1, 1);
    stairs(1:k, y);
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u);
    xlabel('k');
    ylabel('u');
end

%Skok z 0.0 na -0.5 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = -0.5;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    subplot(2, 1, 1);
    stairs(1:k, y);
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u);
    xlabel('k');
    ylabel('u');
end

%Skok z 0.0 na -1 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
u(10:k) = -1;
for i = 7:k
    y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
end

if figures
    subplot(2, 1, 1);
    stairs(1:k, y);
    legend('0,0 na 1,0', '0,0 na 0,5', '0,0 na 0,2', '0,0 na -0,2', '0,0 na -0,5', '0,0 na -1,0');
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u);
    legend('0,0 na 1,0', '0,0 na 0,5', '0,0 na 0,2', '0,0 na -0,2', '0,0 na -0,5', '0,0 na -1,0');
    xlabel('k');
    ylabel('u');
    hold off;
end

if saving
    matlab2tikz('results/2/Answers.tex');
end

%2.2 Charakterystyka statyczna procesu y(u)

% Przygotowanie wektorów dla ch. stat. 
ustat = Umin:0.1:Umax;
ystat = zeros( size(ustat, 1), 1);
iter = 0;

for uk = ustat % dla kolejnych wartości od Umin do Umax symulujemy do stabilizacji i bierzemy y pracy dla danego u pracy
    iter = iter+1;
    y = ones(k, 1)*Ypp; %alokacja wektora o długości symulacji
    u = ones(k,1)*Upp; %Sterowanie stałe równe punktow pracy
    u(10:k) = uk;
    for i = 7:k
        y(i) = symulacja_obiektu7y(u(i-5),u(i-6),y(i-1),y(i-2));
    end
    ystat(iter) = y(k);
end

if figures
    figure;
    plot(ustat, ystat, 'b');
    title('Charakterystyka statyczna');
    xlabel('u');
    ylabel('y');
    legend('Charakterystyka statyczna', 'Location', 'southeast');
end

if saving
    matlab2tikz('results/2/StaticCharacteristic.tex');
end



%% 3. Napisac i omówic program w jezyku MATLAB do symulacji cyfrowego algorytmu PID 
%     oraz algorytmu DMC (w najprostszej wersji analitycznej) dla symulowanego procesu.

%% 4. Dla zaproponowanej trajektorii zmian sygnału zadanego (kilka skoków o róznej wartosci,   
%     przyjac mozliwie duze zmiany punktu pracy, wynikajace z charakterystyki statycznej)
%     dobrac nastawy regulatora PID i parametry algorytmu DMC (dowolna metoda).
%     Omówic metode doboru nastaw i uzasadnic jej zastosowanie. Jakosc regulacji oceniac
%     jakosciowo (na podstawie rysunków przebiegów sygnałów) oraz ilosciowo, wyznaczajac
%     wskaznik jakosci regulacji
%     E = \sum^{k_{konc}}_{k=1} (yzad(k) �?� y(k))^2
%     gdzie k_{konc} oznacza koniec symulacji (zawsze taki sam). Zamiescic wyniki symulacji
%     oraz wartosci wskaznika jakosci E.
mkdir('results/4');
rng(19981998); % ustalenie seedu, żeby rand zawsze dawał to samo przy każdym uruchomieniu
for i = 0:9
    y_zad_traj(1+200*i:200*(i+1)) = Ymin + (Ymax-Ymin)*rand(1);
end
if figures
figure;
plot(1:2000, y_zad_traj, 'r:');
title('Trajektoria zadana');
legend('Trajektoria zadana');
xlabel('k');
ylabel('y');
end

if saving
    matlab2tikz('results/4/TrajektoriaZadana.tex');
end

%% PID

K = 14.36;
Ti = 4.99;
Td = 3.54;

for i = 0:9
    yzad = y_zad_traj(1+200*i);
    if i == 0
        [u, y] = PID(K, Ti, Td, 200, yzad, Upp, Ypp);
    else
        [u, y] = PID(K, Ti, Td, 200, yzad, u_traj(200*i), y_traj(200*i));
    end
    u_traj(1+200*i:200*(i+1)) = u;
    y_traj(1+200*i:200*(i+1)) = y;
end

E = (y_traj-y_zad_traj)*(y_traj-y_zad_traj)';

PIDtitle = sprintf('Algorytm PID K = %g Ti = %g Td = %g E = %g', K, Ti, Td, E);
PIDtitle = strrep(PIDtitle,'.',',');
k = 2000;
if figures
    plotProcess(u_traj, y_traj, PIDtitle);
end
subplot(2,1,1);
hold on;
plot(1:k, y_zad_traj, 'r:');
legend('Wyjście procesu', 'Wartość zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//4//%s.tex', PIDtitle));
end

%% DMC

% D = 140;
% N = 45;
% Nu = 4;
% lambda = 20;
% 
% for i = 0:9
%     yzad = y_zad_traj(1+200*i);
%     if i == 0
%         [u, y] = DMC(s, D, N, Nu, lambda, yzad, 200, Upp,Ypp);
%     else
%         [u, y] = DMC(s, D, N, Nu, lambda, yzad, 200, u_traj(200*i), y_traj(200*i));
%     end
%     u_traj(1+200*i:200*(i+1)) = u;
%     y_traj(1+200*i:200*(i+1)) = y;
% end
% 
% E = (y_traj-y_zad_traj)*(y_traj-y_zad_traj)';
% 
% DMCtitle = sprintf('Algorytm DMC D = %g N = %g Nu = %g lambda = %g  E = %g', D, N, Nu, lambda, E);
% DMCtitle = strrep(DMCtitle,'.',',');
% k = 2000;
% if figures
%     plotProcess(u_traj, y_traj, DMCtitle);
% end
% subplot(2,1,1);
% hold on;
% plot(1:k, y_zad_traj, 'r:');
% legend('Wyj��cie procesu', 'Warto�� zadana', 'Location', 'southwest');
% hold off;
% k = 250;
% 
% if saving
%     matlab2tikz(sprintf('results//4//%s.tex', DMCtitle));
% end

%% 5. W tym samym programie zaimplementowac i omówic rozmyty algorytm PID i rozmyty
%     algorytm DMC w najprostszej wersji analitycznej. Uzasadnic wybór zmiennej,
%     na podstawie której dokonywane jest rozmywanie. Uzasadnic wybór i kształt funkcji
%     przynaleznosci.

% PID 
%zmienna decyzyjna -  do wyboru u(i-1) oraz y(i)
<<<<<<< HEAD
Rvar = 'u(i-1)';
=======
Rvar = 'y(i)';
>>>>>>> origin/Krzychu
[Rmax,Rmin] = setRuleCons(Rvar); 

%funkcja przynale�no�ci - do wyboru gbellmf, gaussmf, trimf, trapmf
memFun = "gbellmf";
nr = 5;  %liczba regulatorow lokalnych/regul
saving = 0;
[Klocal,Tilocal,Tdlocal] = PIDsetLocalParams(nr); 

E = 0;
for i = 0:9
    yzad = y_zad_traj(1+200*i);
    if i == 0
         [u, y] = PID_distributed(Klocal, Tilocal, Tdlocal, 200, yzad, Upp, Ypp, memFun);
    else
        [u, y] = PID_distributed(Klocal, Tilocal, Tdlocal, 200, yzad, u_traj(200*i), y_traj(200*i), memFun);
    end
    u_traj(1+200*i:200*(i+1)) = u;
    y_traj(1+200*i:200*(i+1)) = y;
    
    E = E + sum((yzad * ones(200, 1) - y').^2);
end

PID2title = sprintf('PID dla %d regulator�w lokalnych i funkcji przynale�no�ci %s E = %.2f', nr, memFun, E);
PID2title = strrep(PID2title,'.',',');
savName = sprintf('%d_pid_fun_%s', nr, memFun);
savName = strrep(savName,'.',',');
k = 2000;
if figures
    plotProcess(u_traj, y_traj, PID2title);
end

subplot(2,1,1);
hold on;
plot(1:k, y_zad_traj, 'r:');
legend('Wyj�cie procesu', 'Warto�� zadana', 'Location', 'southwest');
hold off;
k = 250;
if saving
    matlab2tikz(sprintf('results//6//%s.tex', savName));
end



%% DMC 
%zmienna decyzyjna - do wyboru u(i-1) oraz y(i)
<<<<<<< HEAD
Rvar = 'u(i-1)';
%ustawienie zakresu warto�ci zmiennej decyzyjnej
[Rmax,Rmin] = setRuleCons(Rvar); 
%funkcja przynale�no�ci - do wyboru gbellmf, gaussmf, trimf, trapmf
memFun = "gbellmf";

%w funkcji wykonywany jest algorytm rozproszonego regulatora DMC
%argumenty: liczba regulator�w, parametr lambda
testDMC(2,1)
=======
% Rvar = 'u(i-1)';
% [Rmax,Rmin] = setRuleCons(Rvar); 
% %funkcja przynale�no�ci - do wyboru gbellmf, gaussmf, trimf, trapmf
% memFun = "trimf";
% 
% %liczba regulator�w lokalnych 
% nr = 4;
% 
% %zebranie lokalnych odpowiedzi skokowych 
% tau = 10 ; 
% y = ones(k, nr)*Ypp; %alokacja wektora o długości symulacji
% u = ones(k,nr)*Upp; %Sterowanie stałe równe punktow pracy
% yDMC = zeros(k,nr); 
% uDMC = zeros(k,nr); 
% s = zeros(k,nr); 
% 
% centerDistance = (Rmax-Rmin)/(nr-1); 
% centers = Rmin : centerDistance : Rmax; 
% figure
% for reg = 1 : nr
%     u(10:k,reg) = centers(reg);
%     for i = 7:k
%         y(i,reg) = symulacja_obiektu7y(u(i-5,reg),u(i-6,reg),y(i-1,reg),y(i-2,reg));
%     end
%     %przeskalowanie 
%     yDMC(1:end-tau,reg) = y(tau+1:end,reg); 
%     yDMC(end-tau:end,reg) = y(end-tau:end,reg);
%     uDMC(1:end-tau,reg) = u(tau+1:end,reg);
%     uDMC(end-tau:end,reg) = u(end-tau:end,reg);
%     if centers(reg) ~= 0
%         s(1:end,reg) = yDMC(1:end,reg)/centers(reg);
%     end
% %     plot(yDMC(1:end,reg))
% %     hold on
% %     plot(uDMC(1:end,reg))
% %     hold on
% %     figure
%     plot(s(1:end,reg))
%     hold on
%     if centers(reg) ~= 0
%         plot(uDMC(1:end,reg)/centers(reg))
%     end
% end
% 
% 
% D = setD(s);
% N = D;
% Nu = fix(N/3);
% lambda = 1;
% for i = 0:9
%     yzad = y_zad_traj(1+200*i);
%     if i == 0
%         [u, y] = DMC_distributed(s, D, N, Nu, lambda, yzad, 200, Upp,Ypp);
%     else
%         [u, y] = DMC_distributed(s, D, N, Nu, lambda, yzad, 200, u_traj(200*i), y_traj(200*i));
%     end
%     u_traj(1+200*i:200*(i+1)) = u;
%     y_traj(1+200*i:200*(i+1)) = y;
% 
% % 
% %     k = 200 ;
% %     plotProcess(u,y,sprintf('Skok na %.2f',yzad)); 
% %     subplot(2,1,1);
% %     hold on;
% %     plot(1:k, ones(1,k)*yzad, 'r:');
% end
% 
% k = 2000;
% if figures
%     plotProcess(u_traj, y_traj, 'DMC_distributed');
% end
% subplot(2,1,1);
% hold on;
% plot(1:k, y_zad_traj, 'r:');
% legend('Wyj��cie procesu', 'Warto�� zadana', 'Location', 'southwest');
% hold off;
% k = 250;
% 
% 
% 
% 
>>>>>>> origin/Krzychu



%% 6. Dobrac parametry kazdego z lokalnych regulatorów w taki sposób, aby osiagnac mozliwie
%     wysoka jakosc regulacji w okolicach jego punktu pracy (przyjac dla DMC lambda = 1).
%     Wykonac, dla załozonej trajektorii zmian sygnału wartosci zadanej, eksperymenty
%     uwzgledniajac

for nregs = 2 : 10 
    testDMC(nergs,1);
end

%% 7. Dla zaproponowanej trajektorii zmian sygnału zadanego oraz dla róznej liczby regulatorów
%     lokalnych (2, 3, 4, 5, . . . ) spróbowac dobrac parametry lambda dla kazdego z lokalnych
%      regulatorów DMC. Zamiescic wyniki symulacji.


for lambda1  = [5,8,12,16,20]
    for lambda2 = [6,10,14,18,25] 
    testDMC(2,[lambda1,lambda2]);
    end
end

