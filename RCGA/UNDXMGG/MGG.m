function Population = MGG(Param, Population)
% Minimal Generation Gap Selection for UNDX
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p234, 2000 

maxitr = 10;
n_population = Param.n_population;
n_children = Param.n_children;
n_constraint = Param.n_constraint;
n_gene = Param.n_gene;
Pf = Param.Pf;

% ip = randperm(n_population,2);
ip = [ 1 3 ];

c(1,1:n_children) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    for j = 1 : maxitr
%         ip(3) = randi(n_population);
ip(3) = 5;
        if ip(3) ~= ip(1) && ip(3) ~= ip(2)
            break;
        end
    end
    c(i) = getNewChild(Population(ip(1)),Population(ip(2)),Population(ip(3)));
end
c = requestFitnessCalc(Param,c);

f = [ c Population(ip(1:2)) ];

% for i = 1 : n_children+2
%     fprintf(':%e\t%e\t%e\t%e\t%e\t%e\t%e\n', ...
%         i,f(i).gene(1),f(i).gene(2), ...
%         f(i).g(1),f(i).g(2),f(i).f,f(i).phi);
% end

f = SRsort(f,Pf);
% Population(ip(1:2)) = [ f(1),  f(randi([2 n_children+2])) ];
i = floor( rand() * (n_children+1) ) + 2;
Population(ip(1:2)) = [ f(1),  f(i) ];
Population = SRsort(Population,Pf);
