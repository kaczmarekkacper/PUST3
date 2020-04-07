function  plotProcess(u, y, plotTitle)

    global k;
    figure;
    subplot(2, 1, 1);
    stairs(1:k, y, 'b');
    title(plotTitle);
    legend('Wyj�cie procesu');
    xlabel('k');
    ylabel('y');
    subplot(2, 1, 2);
    stairs(1:k, u, 'm');
    legend('Wej�cie procesu');
    %axis([0 2000 0.1 1.5])
    xlabel('k');
    ylabel('u');
end

