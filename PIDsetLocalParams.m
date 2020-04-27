function [Klocal,Tilocal,Tdlocal] = PIDsetLocalParams(nregs)

Klocal = zeros(nregs,1); 
Tilocal = zeros(nregs,1);
Tdlocal = zeros(nregs,1); 

if nregs == 2 
% skok w sterowaniu z -1 na -0.9, Kkryt = 0.35, Tu = 21 * 0.5s
%centrum: -1
Klocal(1) = 0.21 ; Tilocal(1) = 5.25; Tdlocal(1) = 1.31 ; 
% skok w sterowaniu z 0.4 na 0.6, Kkryt = 21.36, Tu = 18 * 0.5s 
%centrum: 1
Klocal(2) = 11 ; Tilocal(2) = 4.5; Tdlocal(2) = 1.125 ; 
end

if nregs == 3
    n = 1 ; %centrum: -1
    Klocal(n) = 0.21 ; Tilocal(n) = 5.25; Tdlocal(n) = 1.31 ; 

	n = 2 ; 
    % centrum -0.35
    % skok w sterowaniu z -0.4 na -0.3,Kkryt = 0.513, Tu = 18 * 0.5s 
    Klocal(n) = 0.3 ; Tilocal(n) = 4.5; Tdlocal(n) = 1.125 ;  
    
    n = 3;
    % centrum 0.5
    %skok w sterowaniu z 0.4 na 0.5, Kkryt = 20.165, Tu = 18 * 0.5s 
    Klocal(n) = 10.5 ; Tilocal(n) = 3; Tdlocal(n) = 0.7 ; 
end

if nregs == 4 
    n = 1 ; %centrum: -1
    Klocal(n) = 0.21 ; Tilocal(n) = 5.25; Tdlocal(n) = 1.31 ; 
    
    n = 2 ; 
    % centrum -0.35
    % skok w sterowaniu z -0.4 na -0.3,Kkryt = 0.513, Tu = 18 * 0.5s 
    Klocal(n) = 0.3 ; Tilocal(n) = 4.5; Tdlocal(n) = 1.125 ;  
    
    n = 3;
    % centrum 0.1
    % skok w sterowaniu z 0.1 na 0.2,Kkryt = 5.02 , Tu =16  * 0.5s 
    Klocal(n) = 3 ; Tilocal(n) = 4; Tdlocal(n) = 1 ;  
    
    n = 4;
    % centrum 0.5
    %skok w sterowaniu z 0.4 na 0.5, Kkryt = 20.165, Tu = 18 * 0.5s 
    Klocal(n) = 10.5 ; Tilocal(n) = 3; Tdlocal(n) = 0.7 ; 
end

if nregs == 5 
    n = 1 ; %centrum: -1
    Klocal(n) = 0.21 ; Tilocal(n) = 5.25; Tdlocal(n) = 1.31 ; 
    
    n = 2 ; 
    % centrum -0.55
    % skok w sterowaniu z -0.6 na -0.5,Kkryt = 0.397, Tu = 18 * 0.5s 
    Klocal(n) = 0.238 ; Tilocal(n) = 4.25; Tdlocal(n) =  1.125;  
    
    n = 3;
    % centrum 0
    % skok w sterowaniu z -0.1 na 0.1,Kkryt =  2.3, Tu = 18 * 0.5s 
    Klocal(n) = 1.38 ; Tilocal(n) = 4.25; Tdlocal(n) = 1.125;  
    
    n = 4;
    % centrum 0.2
    %skok w sterowaniu z 0.1 na 0.3, Kkryt = 6.9, Tu = 18 * 0.5s 
    Klocal(n) =  3; Tilocal(n) = 3; Tdlocal(n) =  1.25; 
    
    n = 5;
    % centrum 0.9
    %skok w sterowaniu z 0.4 na 0.5, Kkryt = 20.165, Tu = 18 * 0.5s 
    Klocal(n) = 10.5 ; Tilocal(n) = 3; Tdlocal(n) = 0.7 ; 
end
%params = [Klocal, Tilocal, Tdlocal]; 

end

