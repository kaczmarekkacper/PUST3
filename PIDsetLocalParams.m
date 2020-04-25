function [Klocal,Tilocal,Tdlocal] = PIDsetLocalParams(nregs)

Klocal = zeros(nregs,1); 
Tilocal = zeros(nregs,1);
Tdlocal = zeros(nregs,1); 

if nregs == 2 
Klocal(1) = 5 ; Tilocal(1) = 0.5; Tdlocal(1) = 7 ; 
Klocal(2) = 9 ; Tilocal(2) = 1; Tdlocal(2) = 2 ; 
end

if nregs == 3
    n = 1 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 0.5; Tdlocal(n) = 7 ; 
    n = 2 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 9 ; Tilocal(n) = 1; Tdlocal(n) = 2 ; 
    n = 3 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 2; Tdlocal(n) = 10 ;
end

if nregs == 4 
    n = 1 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 0.5; Tdlocal(n) = 7 ; 
    n = 2 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 9 ; Tilocal(n) = 1; Tdlocal(n) = 2 ; 
    n = 3 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 2; Tdlocal(n) = 10 ; 
    n = 4 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 10 ; Tilocal(n) = 2; Tdlocal(n) = 5 ; 
    
end

if nregs == 5 
    n = 1 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 0.5; Tdlocal(n) = 7 ; 
    n = 2 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 9 ; Tilocal(n) = 1; Tdlocal(n) = 2 ; 
    n = 3 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 2; Tdlocal(n) = 10 ; 
    n = 4 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 10 ; Tilocal(n) = 2; Tdlocal(n) = 5 ; 
    n = 5 ; %nr regulatora ktorego parametry modyfikujemy
    Klocal(n) = 5 ; Tilocal(n) = 2; Tdlocal(n) = 10 ; 
    
end
%params = [Klocal, Tilocal, Tdlocal]; 

end

