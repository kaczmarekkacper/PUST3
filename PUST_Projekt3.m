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
global k; %ilość próbek symulacji



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
%     E = \sum^{k_{konc}}_{k=1} (yzad(k) − y(k))^2
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

D = 140;
N = 45;
Nu = 4;
lambda = 20;

for i = 0:9
    yzad = y_zad_traj(1+200*i);
    if i == 0
        [u, y] = DMC(s, D, N, Nu, lambda, yzad, 200, Upp,Ypp);
    else
        [u, y] = DMC(s, D, N, Nu, lambda, yzad, 200, u_traj(200*i), y_traj(200*i));
    end
    u_traj(1+200*i:200*(i+1)) = u;
    y_traj(1+200*i:200*(i+1)) = y;
end

E = (y_traj-y_zad_traj)*(y_traj-y_zad_traj)';

DMCtitle = sprintf('Algorytm DMC D = %g N = %g Nu = %g lambda = %g  E = %g', D, N, Nu, lambda, E);
DMCtitle = strrep(DMCtitle,'.',',');
k = 2000;
if figures
    plotProcess(u_traj, y_traj, DMCtitle);
end
subplot(2,1,1);
hold on;
plot(1:k, y_zad_traj, 'r:');
legend('Wyjście procesu', 'Wartość zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//4//%s.tex', DMCtitle));
end

%% 5. W tym samym programie zaimplementowac i omówic rozmyty algorytm PID i rozmyty
%     algorytm DMC w najprostszej wersji analitycznej. Uzasadnic wybór zmiennej,
%     na podstawie której dokonywane jest rozmywanie. Uzasadnic wybór i kształt funkcji
%     przynaleznosci.
%% 6. Dobrac parametry kazdego z lokalnych regulatorów w taki sposób, aby osiagnac mozliwie
%     wysoka jakosc regulacji w okolicach jego punktu pracy (przyjac dla DMC lambda = 1).
%     Wykonac, dla załozonej trajektorii zmian sygnału wartosci zadanej, eksperymenty
%     uwzgledniajac

%% 7. Dla zaproponowanej trajektorii zmian sygnału zadanego oraz dla róznej liczby regulatorów
%     lokalnych (2, 3, 4, 5, . . . ) spróbowac dobrac parametry lambda dla kazdego z lokalnych
%      regulatorów DMC. Zamiescic wyniki symulacji.

