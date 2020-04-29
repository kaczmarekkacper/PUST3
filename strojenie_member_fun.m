name = "gbellmf";
xmin = -1;
xmax = 1;

figure;
funParams = [1,1.5,-1];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
hold on;

funParams = [1,1.5,1];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 

% funParams = [0.15,3.5,0.3];
% mf = fismf(name,funParams); 
% plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
% 
% funParams = [0.4,3.5,0.9];
% mf = fismf(name,funParams); 
% plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
% 
% funParams = [0.4,3.5,0.9];
% mf = fismf(name,funParams); 
% plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
hold off;