function [f, g] = g04(x)


f = 5.3578547 * x(3)^2 + 0.8356891 * x(1) * x(5) + 37.293239 * x(1) - 40792.141;

g(1) =   85.334407 + 0.0056858 * x(2) * x(5) + 0.0006262 * x(1) * x(4) - 0.0022053 * x(3) * x(5) - 92;
g(2) = - 85.334407 - 0.0056858 * x(2) * x(5) - 0.0006262 * x(1) * x(4) + 0.0022053 * x(3) * x(5);
g(3) =   80.51249  + 0.0071317 * x(2) * x(5) + 0.0029955 * x(1) * x(2) + 0.0021813 * x(3) * x(3) - 110;
g(4) = - 80.51249  - 0.0071317 * x(2) * x(5) - 0.0029955 * x(1) * x(2) - 0.0021813 * x(3) * x(3) + 90;
g(5) =   9.300961  + 0.0047026 * x(3) * x(5) + 0.0012547 * x(1) * x(3) + 0.0019085 * x(3) * x(4) - 25;
g(6) = - 9.300961  - 0.0047026 * x(3) * x(5) - 0.0012547 * x(1) * x(3) - 0.0019085 * x(3) * x(4) + 20;