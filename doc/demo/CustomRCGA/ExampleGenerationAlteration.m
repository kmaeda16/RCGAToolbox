function Population = ExampleGenerationAlteration(problem, opts, Population)
% ExampleGenerationAlteration updates population by a simple generation
% alteration algorithm. This function can be used as a template for custom
% generation alteration functions.
% 
% [SYNTAX]
% Population = ExampleGenerationAlteration(problem, opts, Population)
% 
% [INPUT]
% problem    :  Problem structure.
% opts       :  RCGA options. See XXXXXXXXXXX for options.
% Population :  Array of individuals
% 
% [OUTPUT]
% Population :  Array of updated individuals
% 


%% Shortening variable names
n_population = opts.n_population;
n_children = n_population - 1; % Note that opts.n_children is not used.
n_constraint = problem.n_constraint;
n_gene = problem.n_gene;
Pf = opts.Pf;


%% Generating new children
ip = randperm(n_population,2);

c(1,1:n_children) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',zeros,'phi',zeros);

for i = 1 : n_children
    
    % Crossover
    c(i).gene = 0.5 * ( Population(ip(1)).gene + Population(ip(2)).gene );
    
    % Mutation
    mutation_rate = 0.1;
    for j = 1 : n_gene
        if rand < mutation_rate
            c(i).gene(j) = rand;
        end
    end

end
c = RCGArequestFitnessCalc(problem,opts,c);


%% Making a new population
Population = [ Population(1) c ];
Population = RCGAsrsort(Population,Pf);
