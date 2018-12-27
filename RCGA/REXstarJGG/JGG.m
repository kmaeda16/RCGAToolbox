function Population = JGG(problem, opts, Population)

n_population = opts.n_population;
n_parent = opts.n_parent;
Pf = opts.Pf;
selection_type = opts.selection_type;

% Pick up parents from main population
ip = randperm(n_population,n_parent);

p = Population(ip);
c = REXstar(problem,opts,p);

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
