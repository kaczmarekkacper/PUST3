function [Klocal,Tilocal,Tdlocal] = PIDsetLocalParams(nregs)

Klocal = zeros(nregs,1); 
Tilocal = zeros(nregs,1);
Tdlocal = zeros(nregs,1); 

if nregs == 2 
% skok w sterowaniu z -1 na -0.9, Kkryt = 0.35, Tu = 21 * 0.5s
%centrum: -1
Klocal(1) = 0.18 ; Tilocal(1) = 5.5; Tdlocal(1) = 1.45 ; 
% skok w sterowaniu z 0.4 na 0.6, Kkryt = 21.36, Tu = 18 * 0.5s 
%centrum: 1
Klocal(2) = 10 ; Tilocal(2) = 3; Tdlocal(2) = 1.25 ; 
end

if nregs == 3
    n = 1 ; %centrum: -1
    Klocal(n) = 0.18 ; Tilocal(n) = 5.5; Tdlocal(n) = 1.45 ; 

	n = 2 ; 
    % centrum -0.4
    Klocal(n) = 0.25 ; Tilocal(n) = 3; Tdlocal(n) = 1;  
    
    n = 3;
    % centrum 0.9
    Klocal(n) = 10 ; Tilocal(n) = 2.4; Tdlocal(n) = 0.6 ; 
end

if nregs == 4 
    n = 1 ; %centrum: -1
    Klocal(n) = 0.18 ; Tilocal(n) = 5.5; Tdlocal(n) = 1.45 ; 
    
    n = 2 ; 
    % centrum -0.4
    Klocal(n) = 0.25 ; Tilocal(n) = 3; Tdlocal(n) = 1;  
    
    n = 3;
    % centrum 0.3
    Klocal(n) = 3.5 ; Tilocal(n) = 1.7; Tdlocal(n) = 1 ;  
    
    n = 4;
    % centrum 0.9
    Klocal(n) = 10 ; Tilocal(n) = 2.4; Tdlocal(n) = 0.6 ; 
end

if nregs == 5 
    n = 1 ; %centrum: -1
    Klocal(n) = 0.18 ; Tilocal(n) = 5.5; Tdlocal(n) = 1.45 ; 
    
    n = 2 ; 
    % centrum -0.5
    Klocal(n) = 0.21 ; Tilocal(n) = 3.5; Tdlocal(n) =  1.3;  
    
    n = 3;
    % centrum 0
    Klocal(n) = 0.75 ; Tilocal(n) = 2; Tdlocal(n) = 1.2;  
    
    n = 4;
    % centrum 0.2
    Klocal(n) =  3.3; Tilocal(n) = 2.1; Tdlocal(n) =  0.8; 
    
    n = 5;
    % centrum 0.9
    Klocal(n) = 10 ; Tilocal(n) = 2.4; Tdlocal(n) = 0.6 ; 
end
end

