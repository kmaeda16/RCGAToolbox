function Population = JGG(Param, Population)

n_population = Param.n_population;
n_parent = Param.n_parent;
Pf = Param.Pf;
selection_type = Param.selection_type;

% Pick up parents from main population
% ip = randperm(n_population,n_parent);
ip = [1 3 5];

p = Population(ip);

c = REXstar(Param,p);

for j = 1 : Param.n_population
%     fprintf(':%e\t%e\t%e\t%e\t%e\t%e\t%e\n',j,Population(j).gene(1),Population(j).gene(2),Population(j).g(1),Population(j).g(2),Population(j).f,Population(j).phi);
end

switch selection_type
    case 0
        % Chosen from children (Kobayashi, 2009)
        Population(ip) = c(1:n_parent);
    case 1
        % Chosen from family (Kimura et al., 2015)
        f = [p c];
        f = SRsort(f,Pf);
        Population(ip) = f(1:n_parent);
    otherwise
        error('Unexpected selection_type!');
end

Population = SRsort(Population,Pf);
