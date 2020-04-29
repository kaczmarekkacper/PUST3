function [] = drawMF(xmin,xmax,mf)
global funcall ; 
if funcall == 0 
    plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax));
    hold on; 
end
end

