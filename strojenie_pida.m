K = 0.75; % 0.25
Ti = 2; %3
Td = 1.2; %1.125
titleN = 'PID lokalny przy 5 regulatorach lokalnych z centrum w 0';
savName = "3_z_5_pid";
saving = 0;
sim_time = 500;
% ustat, ystat - ch-ka statyczna 21 liczb
yzad = ystat(12);
y_prev = ystat(11);
u_prev = ustat(11);
[u,y] = PID(K,Ti,Td,sim_time, yzad, u_prev, y_prev);
%Rysowanie
global k;
k = sim_time;
E = sum((yzad*ones(sim_time, 1) - y').^2);
titleN = strcat(titleN, sprintf(' E = %.4f', E));
titleN = strrep(titleN,'.',',');
plotProcess(u,y,titleN);
subplot(2, 1, 1);
hold on
stairs(1:sim_time, yzad*ones(sim_time, 1), 'r:');
legend('Wyjœcie procesu', 'Wartoœæ zadana')
hold off;


if saving
    matlab2tikz(sprintf('results//6//%s.tex', savName));
end