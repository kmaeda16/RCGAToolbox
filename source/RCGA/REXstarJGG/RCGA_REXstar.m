function c = RCGA_REXstar(problem, opts, p)
% RCGA_REXstar generates children by using Real-coded Ensemble Crossover
% star (REXstar).
% 
% [SYNTAX]
% c = RCGA_REXstar(problem, opts, p)
% 
% [INPUT]
% problem :  Problem structure:
%            - problem.n_gene: Number of decision variables.
%            - problem.n_constraint: Number of constraint functions. For 
%               unconstained problems, this must be zero.
%            - problem.fitnessfun: Function handle for a fitness function.
%            - problem.decodingfun: Function handle for a decoding 
%               function.
% opts    :  Option structure:
%            - opts.n_population: Population size.
%            - opts.n_children: Number of children.
%            - opts.n_parent: Number of parents.
%            - opts.t_rexstar: Step-size parameter for REXstar/JGG.
%            - opts.selection_type: Selection type for REXstar/JGG 
%               (0 or 1).
%            - opts.Pf: Probability that only the objective function f is 
%               used in comparisons of individuals in the stochastic 
%               ranking.
%            - opts.local: Local optimizer (0 or 1). If it is 1, the local 
%               optimizer is used.
%            - opts.localopts: Options for the local optimizer.
%            - opts.maxgen: Maximum number of generations.
%            - opts.maxtime: Maximum time (sec).
%            - opts.maxeval: Maximum number of fitnessfun evaluations.
%            - opts.vtr: Value to be reached.
%            - opts.n_par: Number of workers in parallel computation.
%            - opts.output_intvl: Interval generation for updating the 
%               transition file and the report file.
%            - opts.out_transition: Name of an output file called the 
%               transition file.
%            - opts.out_best: Name of an output file called the best 
%               individual file.
%            - opts.out_population: Name of an output file called the 
%               final population file.
%            - opts.out_report: Name of an output file called the report 
%               file.
%            - opts.interimreportfun: Function handle for the interim 
%               report function.
%            - opts.finalreportfun: Function handle for the final report 
%               function.
% p       :  Array of parents (chrom structure).
% 
% [OUTPUT]
% c       :  Array of generated children (chrom structure).
% 
% 
% Recommended:
% n_population = 2 * n_gene ~ 20 * n_gene
% n_parent = n_gene + 1
% n_children = 2 * n_gene ~ 3 * n_gene
% t (Step-size parameter) = 2.5 ~ 15
% 
% Required:
% n_parent <= n_children
% 
% See "Kimura S, Sato M, Okada-Hatakeyama M: An Effective Method for the
% Inference of Reduced S-system Models of Genetic Networks. Information and
% Media Technologies 2015, 10(1):166-174.".


%% Shortening variable names
n_parent = opts.n_parent;
n_gene = problem.n_gene;
n_children = opts.n_children;
n_constraint = problem.n_constraint;
t = opts.t_rexstar;
Pf = opts.Pf;


%% Error check
if n_children < n_parent
    error('n_parent <= n_children must be satisfied!');
end


%% Default number of trials of generating children
maxitr = 100; % You can change this line


%% Calculate center of gravity G
G.gene = zeros(1,n_gene);
for i = 1 : n_parent
    G.gene = G.gene + p(i).gene;
end
G.gene = G.gene / n_parent;

%% Initialize children
c(1,1:n_children) = struct('gene',Inf(1,n_gene),'g',Inf(1,n_constraint),'f',Inf,'phi',Inf);


%% Calculate reflection points (which are treated as the children)
count = 0;
index = zeros(1,n_parent);
for i = 1 : n_parent
    c(i).gene = 2.0 * G.gene - p(i).gene;
    if min( 0 <= c(i).gene & c(i).gene <= 1 )
        count = count + 1;
        index(count) = i;
    end
end
index = index(1:count);
c(index) = RCGArequestFitnessCalc(problem,opts,c(index));


%% Store parents and refrection points
r = [p c(1:n_parent)];
r = RCGAsrsort(r,Pf);


%% Calculate center of gravity Gb
Gb.gene = zeros(1,n_gene);
for i = 1 : n_parent
    Gb.gene = Gb.gene + r(i).gene;
end
Gb.gene = Gb.gene / n_parent;


%% Generate children
count = 0;
index = zeros(1,n_children);
for i = 1 : n_children
    flg_eval = 0;
    % If the generated child is out of range, re-generate a child
    for l = 0 : maxitr - 1
        flg_oob = 0;
        if max( c(i).gene < 0 | 1 < c(i).gene )
            flg_oob = 1;
        else
            break;
        end
        % Generation
%         xi = t * rand(1,n_gene);
        for j = 1 : n_gene
            xi(j) = t * rand();
        end
        c(i).gene = G.gene + xi .* ( Gb.gene - G.gene );
        for k = 1 : n_parent
            xi = 2.0 * sqrt( 3.0 / n_parent ) * rand - sqrt( 3.0 / n_parent );
            c(i).gene = c(i).gene + xi .* ( p(k).gene - G.gene );
        end
%         fprintf(':%e\t%e\t%e\t%e\t%e\t%e\t%e\n',i,c(i).gene(1),c(i).gene(2),c(i).g(1),c(i).g(2),c(i).f,c(i).phi);
        flg_eval = 1;
    end
    % If the generated child is out of range, alter the child
    if flg_oob == 1
        c(i).gene( c(i).gene < 0 ) = 0;
        c(i).gene( c(i).gene > 1 ) = 1;
        flg_eval = 1;
    end
    if flg_eval == 1
        count = count + 1;
        index(count) = i;
    end
end

index = index(1:count);
c(index) = RCGArequestFitnessCalc(problem,opts,c(index));


%% Stochastic ranking sort
c = RCGAsrsort(c,Pf);
