function [membership] = MembershipFun(name, xmin,xmax, nregs,x)

%zakres wartoœci, odleg³oœci miêdzy centrami - wykorzystywane przy
%rozk³adzie równomiernym 
% range = xmax - xmin; 
% halfDistance = range/(2*(nregs-1)) ;
% centerDistance = range/(nregs-1); 
% membership = zeros(1,nregs); 
% centers = xmin : centerDistance : xmax; 


global funcall ; %zmienna która warunkuje rysowanie wykresów funkcji przynale¿noœci tylko przy pierwszym wywo³aniu

=======
range = xmax - xmin; 
funParams = []; 
halfDistance = range/(2*(nregs-1)) ;
centerDistance = range/(nregs-1); 
membership = zeros(1,nregs); 
% centers = xmin : centerDistance : xmax; 

>>>>>>> origin/Krzychu
if(nregs == 2)
    centers = [-1;1];
elseif( nregs == 3)
    centers = [-1;-0.4;0.9]; 
elseif(nregs == 4)
    centers = [-1;-0.4;0.3;0.9];
elseif(nregs == 5)
    centers = [-1;-0.55;0;0.2;0.9];
end

if strcmp(name,"gbellmf")
<<<<<<< HEAD
        if funcall == 0
             figure ;
             title(sprintf('Funkcje przynale¿noœci dla %d regulatorów lokalnych',nregs));
             hold on
        end
        if(nregs == 2)
            funParams = [1.4,3.5,-1];
            mf = fismf(name,funParams); 
            drawMF(xmin,xmax,mf); 
            membership(1) = evalmf(mf,x); 
            funParams = [0.5,3.5,1];
            mf = fismf(name,funParams); 
            membership(2) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funcall = 1 ; 
        elseif(nregs == 3)
            funParams = [0.4,3.5,-1];
            mf = fismf(name,funParams); 
            membership(1) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funParams = [0.5,3.5,-0.4];
            mf = fismf(name,funParams); 
            membership(2) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funParams = [0.4,3.5,0.9];
            mf = fismf(name,funParams); 
            membership(3) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funcall = 1;
         elseif(nregs == 4)
            funParams = [0.333,3.5,-1];
            mf = fismf(name,funParams); 
            membership(1) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funParams = [0.4,3.5,-0.4];
            mf = fismf(name,funParams); 
            membership(2) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funParams = [0.15,3.5,0.3];
            mf = fismf(name,funParams); 
            membership(3) = evalmf(mf,x);
            drawMF(xmin,xmax,mf); 
            funParams = [0.4,3.5,0.9];
            mf = fismf(name,funParams); 
            membership(4) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funcall = 1 ;
         elseif(nregs == 5)
            funParams = [0.333,3.5,-1];
            mf = fismf(name,funParams); 
            membership(1) = evalmf(mf,x);
            drawMF(xmin,xmax,mf); 

            funParams = [0.25,3.5,-0.55];
            mf = fismf(name,funParams); 
            membership(2) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 

            funParams = [0.1,3.5,0];
            mf = fismf(name,funParams); 
            membership(3) = evalmf(mf,x);
            drawMF(xmin,xmax,mf); 

            funParams = [0.1,3.5,0.2];
            mf = fismf(name,funParams); 
            membership(4) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 

            funParams = [0.4,3.5,0.9];
            mf = fismf(name,funParams); 
            membership(4) = evalmf(mf,x); 
            drawMF(xmin,xmax,mf); 
            funcall = 1 ;
        end
        if nregs > 5
            for reg  = 1 : nregs
                funParams = [halfDistance,2,centers(reg)];
                mf = fismf(name,funParams); 
                if funcall == 0
                    if reg == 1 
                        figure; 
                    end
                    plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax));
                    hold on

                    if reg == nregs
                        funcall = 1 ; 
                        title(sprintf('Funkcja przynale¿noœci %s dla %d regulatorów lokalnych',name,nregs));
                        legend('off'); 
                        matlab2tikz(sprintf('results//5//DMC//Mf%sRegs%d.tex',name,nregs));
                    end
                end
                membership(reg) = evalmf(mf,x); 
            end
        end
elseif strcmp (name,"gaussmf")
        for reg = 1 : nregs 
            if reg < nregs 
                funParams = [(centers(reg+1)-centers(reg))/2,centers(reg)]; 
            else
                funParams = [(centers(reg)-centers(reg-1))/2,centers(reg)];
            end
            mf = fismf(name,funParams); 
            if funcall == 0
                        if reg == 1 
                            figure; 
                            title(sprintf('Funkcja przynale¿noœci %s dla %d regulatorów lokalnych',name,nregs));
                            legend('off'); 
                            hold on
                        end
                        plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax));
                        hold on
                      
                        if reg == nregs
                            funcall = 1 ; 
                            matlab2tikz(sprintf('results//5//DMC//Memberfun%sRegs%d.tex',name,nregs));
                        end
            end
            membership(reg) = evalmf(mf,x); 
        end
    else

        if strcmp (name,"trimf")
            for reg = 1 : nregs 
                if reg < nregs 
                    if reg > 1 
                        funParams = [centers(reg-1),centers(reg),centers(reg+1)]; 
                    else 
                        funParams = [centers(reg)-1,centers(reg),centers(reg+1)]; 
                    end
                else
                    funParams = [centers(reg-1),centers(reg),centers(reg)+1];
                end
>>>>>>> origin/Krzychu
                mf = fismf(name,funParams); 
                if funcall == 0
                   if reg == 1 
                      figure; 
                      title(sprintf('Funkcja przynale¿noœci %s dla %d regulatorów lokalnych',name,nregs));
                      legend('off'); 
                   end
                   plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax));
                   hold on
                   if reg == nregs
                      funcall = 1 ; 
                      matlab2tikz(sprintf('results//5//DMC//Memberfun%sRegs%d.tex',name,nregs));
                   end
                end
                membership(reg) = evalmf(mf,x); 
            end
        else
            if strcmp (name,"trapmf")
                for reg = 1 : nregs 
                    if reg < nregs 
                        centerDistance2 = centers(reg+1)-centers(reg);
                        if reg > 1 
                            centerDistance1 = centers(reg)-centers(reg-1);
                        else 
                            centerDistance1 = 2;
                        end
                    else
                        centerDistance1 = centers(reg)-centers(reg-1); 
                        centerDistance2 = 2;
                    end
                    
                    funParams = [centers(reg)-centerDistance1*3/4,centers(reg)-centerDistance1/4,centers(reg)+centerDistance2/4,centers(reg)+centerDistance2*3/4]; 
                    mf = fismf(name,funParams); 
                    if funcall == 0
                                if reg == 1 
                                    figure; 
                                    title(sprintf('Funkcja przynale¿noœci %s dla %d regulatorów lokalnych',name,nregs));
                                    legend('off'); 
                                end
                                plot(xmin:0.01:xmax,evalmf(mf,xmin:0.01:xmax));
                                hold on
                                if reg == nregs
                                    funcall = 1 ; 
                                    matlab2tikz(sprintf('results//5//DMC//Mf%sRegs%d.tex',name,nregs));
                                end
                    end
                    membership(reg) = evalmf(mf,x); 
                end
            end
        end  
    end
end

<<<<<<< HEAD


end

=======
>>>>>>> origin/Krzychu

