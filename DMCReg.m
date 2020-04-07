classdef DMCReg < handle
    
    properties
        LAMBDA % macierz lambda
        s % wektor odpowiedzi skokowych
        M % macierz dynamiczna
        Mp % macierz zmian odp skokowej
        K % Wektor wzmocnie�
        D % horyzont dynamiki
        N % horyzont predykcji
        Nu % horyzont sterowania
        u_prev % poprzednie sterowanie
        yzad % warto�� zadana
        deltaup % wektor delta Up
        Umin % minimalna warto�� sterowania
        Umax % maksymalna warto�� sterowniaa
        Umaxchange % max delta sterownaia
    end
    
    methods
        function obj = DMCReg(s, D, N, Nu, lambda, Umin, Umax, Umaxchange)
            % przypisanie do zmiennych obiektu
            obj.N = N;
            obj.D = D;
            obj.Nu = Nu;
            obj.s = s;
            obj.LAMBDA = eye(Nu)*lambda;
            obj.Umin = Umin;
            obj.Umax = Umax;
            obj.Umaxchange = Umaxchange;

            % Wyznaczanie macierzy M
            obj.M = zeros(N, Nu);
            for i=1:1:Nu
                obj.M(i:N,i)=s(1:N-i+1)';
            end

            % Wyznaczanie macierzy Mp
            obj.Mp= zeros(N, D-1);
            for i=1:N
               for j=1:D-1
                  if i+j<=D
                     obj.Mp(i,j)=s(i+j)-s(j);
                  else
                     obj.Mp(i,j)=s(D)-s(j);
                  end   
               end
            end

            % Wektor wzmocnie� - wyznaczany raz (offline)
            obj.K=(obj.M'*obj.M+obj.LAMBDA)\obj.M';
        end
        
        function [] = reset(obj,u_p)
            obj.u_prev = u_p;
            obj.deltaup=zeros(1,obj.D-1)';
        end
        
        function [] = setValue(obj, yzad)
            yzad(1:obj.N)= yzad;
            obj.yzad=yzad';
        end
        
        function u = countValue(obj,y)
            
            % aktualizacja wektora aktualnej warto�ci wyj�cia
            yk=ones(obj.N,1)*y;
            
            % wyliczenie nowego wektora odpowiedzi swobodnej
            y0=yk+obj.Mp*obj.deltaup;
            
            % wyliczenie wektora zmian sterowania
            deltauk=obj.K*(obj.yzad-y0);
            
            %ograniczenie zmian sterowania
            if deltauk(1) > obj.Umaxchange
                deltauk(1) = obj.Umaxchange;
            elseif deltauk(1) < -obj.Umaxchange
                deltauk(1) =-obj.Umaxchange;
            end
            
            % prawo regulacji
            u=obj.u_prev+deltauk(1);
            
            %ograniczenie warto�ci
            if u < obj.Umin
                u = obj.Umin;
            elseif u > obj.Umax
                u = obj.Umax;
            end
            
            %Przepisanie poprzedniej warto�ci
            obj.u_prev = u;
            
            % aktualizacja poprzednich zmian sterowania
            obj.deltaup=[deltauk(1) obj.deltaup(1:end-1)']';
        end
    end
end

