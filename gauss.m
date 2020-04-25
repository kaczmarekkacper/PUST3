function [w] = gauss(x,c,sigma)
w = exp((-(x-c).^2)./(2.*sigma));
end

