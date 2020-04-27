function [membership] = MembershipFun(name, xmin,xmax, nregs,x)
range = xmax - xmin; 
funParams = []; 
halfDistance = range/(2*(nregs-1)) ;
centerDistance = range/(nregs-1); 
membership = zeros(1,nregs); 
centers = xmin : centerDistance : xmax; 
if strcmp(name,"gbellmf")
    if(nregs == 2)
        funParams = [1.4,3.5,-1];
        mf = fismf(name,funParams); 
        membership(1) = evalmf(mf,x); 
        funParams = [0.5,3.5,1];
        mf = fismf(name,funParams); 
        membership(2) = evalmf(mf,x); 
    elseif(nregs == 3)
        funParams = [0.4,3.5,-1];
        mf = fismf(name,funParams); 
        membership(1) = evalmf(mf,x); 
        funParams = [0.5,3.5,-0.4];
        mf = fismf(name,funParams); 
        membership(2) = evalmf(mf,x); 
        funParams = [0.4,3.5,0.9];
        mf = fismf(name,funParams); 
        membership(3) = evalmf(mf,x); 
    elseif(nregs == 4)
        funParams = [0.333,3.5,-1];
        mf = fismf(name,funParams); 
        membership(1) = evalmf(mf,x); 
        funParams = [0.4,3.5,-0.4];
        mf = fismf(name,funParams); 
        membership(2) = evalmf(mf,x); 
        funParams = [0.15,3.5,0.3];
        mf = fismf(name,funParams); 
        membership(3) = evalmf(mf,x); 
        funParams = [0.4,3.5,0.9];
        mf = fismf(name,funParams); 
        membership(4) = evalmf(mf,x); 
    elseif(nregs == 5)
        funParams = [0.333,3.5,-1];
        mf = fismf(name,funParams); 
        membership(1) = evalmf(mf,x); 
        
        funParams = [0.25,3.5,-0.55];
        mf = fismf(name,funParams); 
        membership(2) = evalmf(mf,x); 
        
        funParams = [0.1,3.5,0];
        mf = fismf(name,funParams); 
        membership(3) = evalmf(mf,x); 
        
        funParams = [0.1,3.5,0.2];
        mf = fismf(name,funParams); 
        membership(4) = evalmf(mf,x); 
            
        funParams = [0.4,3.5,0.9];
        mf = fismf(name,funParams); 
        membership(4) = evalmf(mf,x); 
    else
        for reg  = 1 : nregs
            funParams = [halfDistance,1.5,centers(reg)];
            mf = fismf(name,funParams); 
            membership(reg) = evalmf(mf,x); 
        end
    end
elseif strcmp (name,"gaussmf")
    for reg = 1 : nregs 
        funParams = [range/(5*centerDistance),centers(reg)]; 
        mf = fismf(name,funParams); 
%             plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%             hold on
        membership(reg) = evalmf(mf,x); 
    end
elseif strcmp (name,"trimf")
    for reg = 1 : nregs 
        funParams = [centers(reg)-centerDistance,centers(reg),centers(reg)+centerDistance]; 
        mf = fismf(name,funParams); 
%                 plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%                 hold on
        membership(reg) = evalmf(mf,x); 
    end
elseif strcmp (name,"trapmf")
    for reg = 1 : nregs 
        funParams = [centers(reg)-centerDistance*3/4,centers(reg)-centerDistance/4,centers(reg)+centerDistance/4,centers(reg)+centerDistance*3/4]; 
        mf = fismf(name,funParams); 
%                     plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%                     hold on
        membership(reg) = evalmf(mf,x); 
    end
end
end
