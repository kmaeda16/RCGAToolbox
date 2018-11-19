function c = UNDX(p1, p2, p3)
% Function of "Unimordal Distribution Crossover"
% p1 and p2 are main parents, and p3 is sub parent.
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p261, 2000

a = 0.5;
b = 0.35;

if length(p1) ~= length(p2) || length(p2) ~= length(p3) || length(p3) ~= length(p1)
    error('Lengths of p1, p2 and p3 must be the same!');
end

n = length(p1.gene);
d1 = norm(p2.gene-p1.gene);
d2 = norm((p3.gene-p1.gene)-(p3.gene-p1.gene)*(p2.gene-p1.gene)'/((p2.gene-p1.gene)*(p2.gene-p1.gene)')*(p2.gene-p1.gene));
e1 = (p2.gene-p1.gene)/norm(p2.gene-p1.gene);

t = randn(1,n)*b*d2/n^0.5;
t = t - (t*e1')*e1;
t = t + randn*a*d1*e1;
c.gene = (p1.gene+p2.gene)/2 + t;
