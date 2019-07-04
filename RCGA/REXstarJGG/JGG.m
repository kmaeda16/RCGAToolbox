function Population = JGG(problem, opts, Population)
% JGG updates population by using Just Generation Gap
% 
% [SYNTAX]
% Population = JGG(problem, opts, Population)
% 
% [INPUT]
% problem   :  Problem structure.
% opts      :  RCGA options. See XXXXXXXXXXX for options.
% Population:  Array of individuals
% 
% [OUTPUT]
% Population:  Array of updated individuals


%% Shortening variable names
n_population = opts.n_population;
n_parent = opts.n_parent;
Pf = opts.Pf;
selection_type = opts.selection_type;


%% Pick up parents from main population
ip = randperm(n_population,n_parent);


%% Generating children
p = Population(ip);
c = REXstar(problem,opts,p);


%% Updating population
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


%% Stochastic ranking sort
Population = SRsort(Population,Pf);
