function [ruleMax,ruleMin] = setRuleCons(ruleVariable)
global Umax 
global Umin
global Ymax
global Ymin
if strcmp(ruleVariable,'u(i-1)') 
    ruleMax = Umax; 
    ruleMin = Umin;
elseif strcmp(ruleVariable,'y(i)')
    ruleMax = Ymax;
    ruleMin = Ymin; 
end


end

