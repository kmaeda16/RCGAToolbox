function c = UNDX(p1, p2, p3)
% Function of "Unimordal Distribution Crossover"
% p1 and p2 are main parents, and p3 is sub parent.
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p261, 2000

alpha = 0.5;
beta = 0.35;

if length(p1) ~= length(p2) || length(p2) ~= length(p3) || length(p3) ~= length(p1)
    error('Lengths of p1, p2 and p3 must be the same!');
end

% %{
n = length(p1.gene);

% v12 is a vector from p1 to p2, v13 is a vector from p1 to p3
v12 = p2.gene - p1.gene;
v13 = p3.gene - p1.gene;

% temp1 is an inner product (v12,v13), temp2 is (v12,v12)
temp1 = dot(v12,v13);
temp2 = dot(v12,v12);
d1 = sqrt(temp2); % =  norm(v12);

if temp2 == 0 || d1 == 0
    d2 = norm(v13);
    e1 = zeros(1,n);
else
    d2 = norm( temp1 / temp2 * v12 - v13 );
    e1 = v12 / d1;
end

t = randn(1,n) * beta * d2 / sqrt(n);
% t = ones(1,n) * beta * d2 / sqrt(n);
temp1 = dot(t,e1);
temp2 = alpha * d1 * randn;
% temp2 = alpha * d1 * 1;
c.gene = 0.5 * ( p1.gene + p2.gene ) + t + ( temp2 - temp1 ) * e1;

% t = randn(1,n) * beta * d2 / sqrt(n);
% temp1 = dot(t,e1);
% temp2 = alpha * d1 * randn;
% c.gene = 0.5 * ( p1.gene + p2.gene ) + t + ( alpha * d1 * randn - dot(t,e1) ) * e1;

%}

%{
n = length(p1.gene);
d1 = norm( p2.gene - p1.gene );
d2 = norm( (p3.gene - p1.gene) - (p3.gene - p1.gene) * (p2.gene - p1.gene)' ...
    / ( (p2.gene - p1.gene) * (p2.gene - p1.gene)' ) * (p2.gene - p1.gene) );
e1 = (p2.gene-p1.gene) / norm(p2.gene-p1.gene);

t = randn(1,n) * beta * d2 / n^0.5;
% t = ones(1,n) * beta * d2 / n^0.5;
t = t - ( t * e1' ) * e1;
t = t + randn * alpha * d1 * e1;
% t = t + 1 * alpha * d1 * e1;
c.gene = ( p1.gene + p2.gene ) / 2 + t;


% c.gene = ( p1.gene + p2.gene ) / 2 + t - ( t * e1' ) * e1 + randn * alpha * d1 * e1;
%}

for i = 1 : length(c.gene)
    if isnan(c.gene(i))
        fprintf('break\n');
    end
end
