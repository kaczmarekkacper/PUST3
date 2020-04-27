name = "gbellmf";
xmin = -1;
xmax = 1;

figure;
funParams = [0.333,3.5,-1];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
hold on;

funParams = [0.25,3.5,-0.55];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 

funParams = [0.1,3.5,0];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 

funParams = [0.1,3.5,0.2];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 

funParams = [0.4,3.5,0.9];
mf = fismf(name,funParams); 
plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
hold off;