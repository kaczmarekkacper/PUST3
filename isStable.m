function [stable] = isStable(data,epsilon)

diff1 = max(data) - mean(data); 
diff2 = mean(data) - min(data);
if diff1 > epsilon ||  diff2 > epsilon
    stable = 0 ;
else
    stable = 1 ;
end


end

