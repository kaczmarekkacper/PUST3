function [D] = setD(s)
D = zeros(length(s(1,1:end)),1); 
eps = 0.00001 ; 
for i  = 1 : length(s(1,1:end))
    for j  = 1 : length (s(1:end,1))
        if isStable(s(j:end,i),eps)
           D(i) = j + 20 ;
           break;
        end
    end
end
end

