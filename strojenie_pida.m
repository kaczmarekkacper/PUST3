K = 1.38; % 4.14
Ti = 4; %4.25
Td = 1.125; %1.125
sim_time = 500;
% ustat, ystat - ch-ka statyczna 21 liczb
yzad = ystat(12);
y_prev = ystat(10);
u_prev = ustat(10);
[u,y] = PID(K,Ti,Td,sim_time, yzad, u_prev, y_prev);
%Rysowanie
global k;
k = sim_time;
plotProcess(u,y,'');
subplot(2, 1, 1);
hold on
stairs(1:sim_time, yzad*ones(sim_time, 1), 'm');
hold off;

hLeg = legend('example');
set(hLeg,'visible','off');

saving = 0;
if saving
    matlab2tikz(sprintf('results//6//%s.tex', 'member_f_3'));
end