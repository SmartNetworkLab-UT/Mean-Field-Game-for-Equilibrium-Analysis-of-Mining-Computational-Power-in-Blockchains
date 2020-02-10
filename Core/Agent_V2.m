classdef Agent_V2 < handle

    properties
        x_min;
        x_max;
        A;
        p;
        gamma;
        initial_gamma;
        adaptingGamma
        x;
        x_f_BR;
        u_f_BR;
        x_f_MF;
        u_f_MF;
        initial_x;
        
    end

    methods 

        function obj = Agent_V2(x_min, x_max, A, p, gamma, x)
            obj.x_min = x_min;
            obj.x_max = x_max;
            obj.A = A;
            obj.p = p;
            obj.gamma = gamma;
            obj.initial_gamma = gamma;
            obj.adaptingGamma = 1;
            obj.x = x;
            obj.initial_x = x;
        end
        
        function RollBack(obj)
           obj.x = obj.initial_x; 
           obj.gamma = obj.initial_gamma;
        end
        
        function DiffX_BR = UpdateStrategyBR(obj, Sum_All)
            initX = obj.x;
            S_i = Sum_All - obj.x;
            obj.gamma = obj.gamma*obj.adaptingGamma;
            X_BR = (obj.A*(S_i/obj.p))^0.5 - S_i;
            if X_BR<obj.x_min
                X_BR = obj.x_min;
            elseif X_BR>obj.x_max
                X_BR = obj.x_max;
            end
            obj.x = (1 - obj.gamma)*obj.x + obj.gamma*X_BR;
            X_BR = obj.x;
            DiffX_BR = X_BR - initX;
        end
        
        function DiffX_MF = UpdateStrategyMF(obj, Sum_All)
            initX = obj.x;
            obj.x = obj.x + obj.gamma*(obj.A/Sum_All - obj.p);
            X_MF = obj.x;
            obj.gamma = obj.gamma*obj.adaptingGamma;
            
            if X_MF<obj.x_min
                X_MF = obj.x_min;
            elseif X_MF>obj.x_max
                X_MF = obj.x_max;
            end
            obj.x = X_MF;
            DiffX_MF = X_MF - initX;
        end
        
        function Util = UtilityCalc(obj, Sum_All)
           Util = obj.A*(obj.x/Sum_All) - obj.p*obj.x; 
        end
    end    
end

