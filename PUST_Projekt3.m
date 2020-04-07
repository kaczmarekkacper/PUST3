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
global k; %iloœæ próbek symulacji



Upp = 0;
Ypp = 0;
Umin = -1;
Umax = 1;
Ymax = 2.721; % wyznaczone eksperymentalnie
Ymin = 1.279;
T = 0.5;
k = 250;

figures = 1; % czy wyœwietlaæ wykresy
saving = 1; % czy zapisywaæ dane

%% 1. Sprawdzic poprawnosc wartosci Upp, Ypp.

mkdir('results/1');

y = zeros(k, 1); %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
for i = 3:k
    y(i) = symulacja_obiektu7y(Upp,Upp,y(i-1),y(i-2));
end

%rysowanie wykresów
if figures
    plotProcess(u, y,'Sprawdzenie punktu pracy');
end

%zapisanie plików do tex, aby wrzuciæ wykres w tikz do latex
if saving
    matlab2tikz('results/1/CheckWorkpoint.tex');
end

%% 2. Wyznaczyc symulacyjnie odpowiedzi skokowe procesu dla kilku zmian sygna³u sterujacego,
%     przy uwzglednieniu ograniczen wartosci tego sygna³u, jego wartosc na poczatku
%     eksperymentu wynosi Upp. Narysowac te odpowiedzi na jednym rysunku. Narysowac
%     charakterystyke statyczna procesu y(u). Czy w³asciwosci statyczne i dynamiczne procesu
%     sa (w przyblizeniu) liniowe? Jezeli tak, okreslic wzmocnienie statyczne procesu.

mkdir('results/2');

% 2.1 Odpowiedzi skokowe

%Skok z 0 na 1 w chwili 10
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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
y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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

for uk = ustat % dla kolejnych wartoœci od Umin do Umax symulujemy do stabilizacji i bierzemy y pracy dla danego u pracy
    iter = iter+1;
    y = ones(k, 1)*Ypp; %alokacja wektora o d³ugoœci symulacji
    u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
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



%% 3. Przekszta³cic jedna z otrzymanych odpowiedzi w taki sposób, aby otrzymac odpowiedz
%     skokowa wykorzystywana w algorytmie DMC, tzn. zestaw liczb s1, s2, . . . (przy skoku
%     jednostkowym sygna³u sterujacego: od chwili k = 0 w³acznie sygna³ sterujacy ma
%     wartosc 1, w przesz³osci jest zerowy). Zamiescic rysunek odpowiedzi skokowej.

tau = 10; %opoznienie

%przygotowanie wektorów
yDMC = ones(k, 1); %alokacja wektora o d³ugoœci symulacji
uDMC = ones(k,1); %Sterowanie sta³e równe punktow pracy

%Skok z 0.8 na 1.5 w chwili 10
y = ones(k, 1)*2; %alokacja wektora o d³ugoœci symulacji
u = ones(k,1)*Upp; %Sterowanie sta³e równe punktow pracy
u(10:k) = Umax;
for i = 12:k
    y(i) = symulacja_obiektu7Y(u(i-10),u(i-11),y(i-1),y(i-2));
end

if figures
    plotProcess(u, y, 'Model DMC przed znormalizowaniem');
end

if saving
    matlab2tikz('results/3/DMCBefore.tex');
end

%sprowadzenie odp skokowej do normy
yDMC(1:end-tau) = y(tau+1:end); yDMC(end-tau:end) = y(end-tau:end);
uDMC(1:end-tau) = u(tau+1:end); uDMC(end-tau:end) = u(end-tau:end);
yDMC = (yDMC-2)/(1.5-0.8);
uDMC = (uDMC-0.8)/(1.5-0.8);

%wektor odpowiedzi skokowych
global s;
s = yDMC;

if figures
    plotProcess(uDMC, yDMC, 'Model DMC po znormalizowaniu');
end

if saving
    matlab2tikz('results/3/DMCAfter.tex');
end

%% 4. Napisac i omówic program w jezyku Matlab do symulacji cyfrowego algorytmu PID
%     oraz algorytmu DMC (w najprostszej wersji analitycznej) dla symulowanego procesu.
%     Istniejace ograniczenia wartosci sygna³u sterujacego oraz ograniczenie szybkosci zmian
%     tego sygna³u
%     -/\Umax <= /\U(k) <= /\Umax
%     gdzie /\Umax = 0,2, uwzglednic odpowiednio ograniczajac (przycinajac) wyznaczony
%     przez regulator sygna³ sterujacy.
mkdir('results/4');
yzad = Ypp+0.1;

%% 4.1 PID

K = 3.583;
Ti = 1000000;
Td = 0;
[u,y] = PID(K, Ti, Td, k, yzad, Upp, Ypp);

%PIDtitle = sprintf('Algorytm PID K = %g Ti = %g Td = %g', K, Ti, Td);
%PIDtitle = sprintf('Szukanie wzmocnienia statycznego K = %g', K);
PIDtitle = 'Zbli¿enie na oscylacje';
PIDtitle = strrep(PIDtitle,'.',',');
if figures
    plotProcess(u, y, PIDtitle);
end
subplot(2,1,1);
hold on;
plot(1:k, ones(1,k)*yzad, 'r:');
legend('Wyjœcie procesu', 'Wartoœæ zadana', 'Location', 'southwest');
%axis([0 2000 1.7 2.2])
hold off;
if saving
    matlab2tikz(sprintf('results//4//%s.tex', PIDtitle));
end
%% 4.2 DMC


%Paramtery domyœlne
D = 250;
N = D;
Nu = N;
lambda = 20;

% symulacja DMC
[u, y] = DMC(s, D, N, Nu, lambda, yzad, k, Upp, Ypp);

DMCtitle = sprintf('Algorytm DMC D = %g N = %g Nu = %g lambda = %g', D, N, Nu, lambda);
DMCtitle = strrep(DMCtitle,'.',',');
if figures
    plotProcess(u, y, DMCtitle);
end
subplot(2,1,1);
hold on;
plot(1:k, ones(1,k)*yzad, 'r:');
legend('Wyjœcie procesu', 'Wartoœæ zadana');
hold off;

if saving
    matlab2tikz(sprintf('results//4//%s.tex', DMCtitle));
end

%% 5. Dla zaproponowanej trajektorii zmian sygna³u zadanego (kilka skoków o róznej amplitudzie)
%     dobrac nastawy regulatora PID i parametry algorytmu DMC metoda eksperymentalna.
%     Jakosc regulacji oceniac jakosciowo (na podstawie rysunków przebiegów
%     sygna³ów) oraz ilosciowo, wyznaczajac wskaznik jakosci regulacji
%     E = \sum^{k_{konc}}_{k=1} (yzad(k) ? y(k))2
%     gdzie kkonc oznacza koniec symulacji (zawsze taki sam). Zamiescic wybrane wyniki
%     symulacji (przebiegi sygna³ów wejsciowych i wyjsciowych procesu oraz wartosci wskaznika E).

mkdir('results/5');

%% Trajektoria zadana

u_traj = zeros(1, 2000);
y_traj = zeros(1, 2000);
y_zad_traj = zeros(1, 2000);

rng(19981998); % ustalenie seedu, ¿eby rand zawsze dawa³ to samo przy ka¿dym uruchomieniu

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
    matlab2tikz('results/5/TrajektoriaZadana.tex');
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
legend('Wyjœcie procesu', 'Wartoœæ zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//5//%s.tex', PIDtitle));
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
legend('Wyjœcie procesu', 'Wartoœæ zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//5//%s.tex', DMCtitle));
end

%% 6. Dla zaproponowanej trajektorii zmian sygna³u zadanego dobrac nastawy regulatora
%     PID i parametry algorytmu DMC (N, Nu, lambda) w wyniku optymalizacji wskaznika jakosci
%     regulacji E. Omówic dobór parametrów optymalizacji. Zamiescic wyniki symulacji dla
%     optymalnych regulatorów.
mkdir('results/6');

%% PID
 
K0 = 13;
Ti0 = 3;
Td0 = 3;
x0 = [K0,Ti0,Td0]; 
A = [-1 0 0;
    0 -1 0;
    0 0 -1];
b = [0;
    0;
    0];
PIDOptimalParameters = fmincon(@PIDFunctionToMinimize, x0, A, b); 
K = PIDOptimalParameters(1); 
Ti = PIDOptimalParameters(2); 
Td = PIDOptimalParameters(3);

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
legend('Wyjœcie procesu', 'Wartoœæ zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//6//%s.tex', PIDtitle));
end


%% DMC

A = [-1 1 0 0;
    0 -1 1 0];
b = [0;
    0];
DMCOptimalParameters = ga(@DMCFunctionToMinimize, 4, A, b, [], [], [140, 1, 1, 10e-5], [250, 250, 250, 1000], [], 1:3);
D = DMCOptimalParameters(1); 
N = DMCOptimalParameters(2); 
Nu = DMCOptimalParameters(3); 
lambda = DMCOptimalParameters(4);

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
legend('Wyjœcie procesu', 'Wartoœæ zadana', 'Location', 'southwest');
hold off;
k = 250;

if saving
    matlab2tikz(sprintf('results//6//%s.tex', DMCtitle));
end

