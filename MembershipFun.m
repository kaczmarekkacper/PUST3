function [membership] = MembershipFun(name, xmin,xmax, nregs,x)
range = xmax - xmin; 
funParams = []; 
halfDistance = range/(2*(nregs-1)) ;
centerDistance = range/(nregs-1); 
membership = zeros(1,nregs); 
centers = xmin : centerDistance : xmax; 
if strcmp(name,"gbellmf")
    for reg  = 1 : nregs
        funParams = [halfDistance,1.5,centers(reg)];
        mf = fismf(name,funParams); 
%         if x == 0
%             plot(xmin:0.1:xmax,evalmf(mf,xmin:0.1:xmax)); 
%             hold on
%         end
        membership(reg) = evalmf(mf,x); 
    end
else 

    if strcmp (name,"gaussmf")
        for reg = 1 : nregs 
            funParams = [range/(5*centerDistance),centers(reg)]; 
            mf = fismf(name,funParams); 
%             plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%             hold on
            membership(reg) = evalmf(mf,x); 
        end
    else

        if strcmp (name,"trimf")
            for reg = 1 : nregs 
                funParams = [centers(reg)-centerDistance,centers(reg),centers(reg)+centerDistance]; 
                mf = fismf(name,funParams); 
%                 plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%                 hold on
                membership(reg) = evalmf(mf,x); 
            end
        else

            if strcmp (name,"trapmf")
                for reg = 1 : nregs 
                    funParams = [centers(reg)-centerDistance*3/4,centers(reg)-centerDistance/4,centers(reg)+centerDistance/4,centers(reg)+centerDistance*3/4]; 
                    mf = fismf(name,funParams); 
%                     plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax)); 
%                     hold on
                    membership(reg) = evalmf(mf,x); 
                end
            end
        end
        
    end
end

end


