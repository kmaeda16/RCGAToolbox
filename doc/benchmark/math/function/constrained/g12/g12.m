function [f, g] = g12(x)
% Constrained benchmark function g12.
% 
% The optimum is located at x* = (5, 5, 5) where f(x*) = -1.
% 
% [SYNTAX]
% [f, g] = g12(x)
% 
% [INPUT]
% x : Decision variables (3 dimensional)
% 
% [OUTPUT]
% f : Objective function
% g : Constraint functions (1 dimensional)


f = - ( 100.0 - ( x(1) - 5 ) ^ 2 - ( x(2) - 5 ) ^ 2 - ( x(3) - 5 ) ^ 2  ) / 100.0;

flg = 0;
for p = 1 : 9
    for q = 1 : 9
        for r = 1 : 9
            if flg == 0
                csn = ( x(1) - p ) ^ 2 + ( x(2) - q ) ^ 2 + ( x(3) - r ) ^ 2 - 0.0625;
                if csn <= 0
                    g(1) = csn;
                    flg = 1;
                end
            end
        end
    end
end

if flg == 0
    g(1) = 1;
end
