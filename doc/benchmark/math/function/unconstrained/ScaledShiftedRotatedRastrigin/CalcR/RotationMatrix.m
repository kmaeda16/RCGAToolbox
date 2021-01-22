function [f, g] = RotationMatrix(x)

global random_vs;

n = length(x);
R = reshape(x,sqrt(n),sqrt(n));

f = 0;
for i = 1 : length(random_vs);
    v = random_vs(:,i);
    f = f + ( norm(R * v) - norm(v) ) ^ 2;
end

g = -det(R);
