function Population = MGG(problem, opts, Population)
% MGG updates population by using Minimal Generation Gap
% 
% [SYNTAX]
% Population = MGG(problem, opts, Population)
% 
% [INPUT]
% problem   :  Problem structure.
% opts      :  RCGA options. See XXXXXXXXXXX for options.
% Population:  Array of individuals
% 
% [OUTPUT]
% Population:  Array of updated individuals
% 
% 
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p234, 2000 


%% Shortening variable names
n_population = opts.n_population;
n_children = opts.n_children;
n_constraint = problem.n_constraint;
n_gene = problem.n_gene;
Pf = opts.Pf;


%% Default number of trials of generating children
maxitr = 10; % You can change this line


%% Generating new children
ip = randperm(n_population,2);

c(1,1:n_children) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    for j = 1 : maxitr
        ip(3) = randi(n_population);
        if ip(3) ~= ip(1) && ip(3) ~= ip(2)
            break;
        end
    end
    c(i) = getNewChild(Population(ip(1)),Population(ip(2)),Population(ip(3)));
end
c = requestFitnessCalc(problem,opts,c);


%% Making a family
f = [ c Population(ip(1:2)) ];
f = SRsort(f,Pf);


%% Making a new population
Population(ip(1:2)) = [ f(1),  f(randi([2 n_children+2])) ];
Population = SRsort(Population,Pf);
