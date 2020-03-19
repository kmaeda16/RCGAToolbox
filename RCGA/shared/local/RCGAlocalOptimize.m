function [improvedChrom, localneval] = RCGAlocalOptimize(problem, opts, chrom)
% RCGAlocalOptimize improves chrom by a local optimizer fmincon. Please
% note that Optiimzation Toolbox is required to use fmincon and thus
% RCGAlocalOptimize.
% 
% [SYNTAX]
% [improvedChrom, localneval] = RCGAlocalOptimize(problem, opts, chrom)
% 
% [INPUT]
% problem       :  Problem structure
% opts          :  RCGA options. See XXXXXXXXXXX for options.
% chrom         :  Individual
% 
% [OUTPUT]
% improvedChrom :  Improved individual


localoptimopts = opts.localoptimopts;

LB = zeros(1,problem.n_gene); % Lower bound
UB = ones(1,problem.n_gene);  % Upper bound
    
ObjectiveFunction  = @(gene) obj_wrapper(problem.fitnessfun, problem.decodingfun, gene);
ConstraintFunction = @(gene) cst_wrapper(problem.fitnessfun, problem.decodingfun, gene);

if problem.n_constraint == 0
    [improvedChrom.gene, ~, ~, output] = fmincon(ObjectiveFunction,chrom.gene,[],[],[],[],LB,UB,[],localoptimopts);
else
    [improvedChrom.gene, ~, ~, output] = fmincon(ObjectiveFunction,chrom.gene,[],[],[],[],LB,UB,ConstraintFunction,localoptimopts);
end

[ improvedChrom.f, improvedChrom.g, improvedChrom.phi ] = RCGAgetFitness(problem,improvedChrom);

% If improvedChrom is worse than chrom, then chrom is returned.
if ~( improvedChrom.phi < chrom.phi || ( improvedChrom.phi == chrom.phi && improvedChrom.f < chrom.f ) )
    improvedChrom = chrom;
end

localneval = output.funcCount + 1;