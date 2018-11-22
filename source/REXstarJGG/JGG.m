function Population = JGG(Param, Population)

n_population = Param.n_population;
n_parent = Param.n_parent;
Pf = Param.Pf;
selection_type = Param.selection_type;

% Pick up parents from main population
ip = randperm(n_population,n_parent);

p = Population(ip);

c = REXstar(Param,p);

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
