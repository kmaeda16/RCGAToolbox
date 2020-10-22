function c = RCGA_UNDX(p1, p2, p3)
% RCGA_UNDX genertes a child by using Unimordal Distribution Crossover
% (UNDX).
% 
% [SYNTAX]
% c = RCGA_UNDX(p1, p2, p3)
% 
% [INPUT]
% p1 :  First main parent (chrom structure).
% p2 :  Second main parent (chrom structure).
% p3 :  Sub-parent (chrom structure).
% 
% [OUTPUT]
% c  :  Generated child (chrom structure).
% 
% 
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p261, 2000


%% Error check
if length(p1.gene) ~= length(p2.gene) || ...
        length(p2.gene) ~= length(p3.gene) || ...
        length(p3.gene) ~= length(p1.gene)
    error('Lengths of p1, p2 and p3 must be the same!');
end


%% These are recommened values
alpha = 0.5;
beta = 0.35;


%% Getting the number of genes
n = length(p1.gene);


%% Generating a child
v12 = p2.gene - p1.gene; % Vector from p1 to p2
v13 = p3.gene - p1.gene; % Vector from p1 to p3
temp1 = dot(v12,v13); % Inner product (v12,v13)
temp2 = dot(v12,v12); % Inner product (v12,v12)
d1 = sqrt(temp2); % =  norm(v12);

if temp2 == 0 || d1 == 0
    d2 = norm(v13);
    e1 = zeros(1,n);
else
    d2 = norm( temp1 / temp2 * v12 - v13 );
    e1 = v12 / d1;
end

% t = randn(1,n) * beta * d2 / sqrt(n);
for i = 1 : n
    t(i) = randn() * beta * d2 / sqrt(n);
end
temp1 = dot(t,e1);
temp2 = alpha * d1 * randn();
% temp2 = alpha * d1 * randn_test();
c.gene = 0.5 * ( p1.gene + p2.gene ) + t + ( temp2 - temp1 ) * e1;
