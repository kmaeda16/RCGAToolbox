function Population = MGG(Param, Population)
% Minimal Generation Gap Selection for UNDX
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p234, 2000 

maxitr = 10;
n_population = Param.n_population;
n_children = Param.n_children;
n_constraint = Param.n_constraint;
n_gene = Param.n_gene;
Pf = Param.Pf;

ip = randperm(n_population,2);

f(1,1:n_children+2) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    for j = 1 : maxitr
        ip(3) = randi(n_population);
        if ip(3) ~= ip(1) && ip(3) ~= ip(2)
            break;
        end
    end
    f(i) = getNewChild(Population(ip(1)),Population(ip(2)),Population(ip(3)));
end
f = requestFitnessCalc(Param,f);

f = [ f, Population(ip(1:2)) ];
f = SRsort(f,Pf);
Population(ip(1:2)) = [ f(1),  f(randi([2 n_children+2])) ];
Population = SRsort(Population,Pf);
