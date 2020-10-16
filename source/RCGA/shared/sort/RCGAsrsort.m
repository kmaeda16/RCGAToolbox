function Sorted_Population = RCGAsrsort(Population, Pf)
% SRsort is a function of stochastic ranking sort.
% 
% [SYNTAX]
% Population = RCGAsrsort(Population, Pf)
% 
% [INPUT]
% Population :  (Unsorted) Population
% Pf         :  Probability that only f is used in comparisons of
%               individuals
% 
% [OUTPUT]
% Population :  Sorted population


%% Preparation
n_population = length(Population);
n_gene = length(Population(1).gene);
n_constraint = length(Population(1).g);

index_tmp = zeros(1,n_population);
f_tmp = zeros(1,n_population);
phi_tmp = zeros(1,n_population);

for i = 1 : n_population
    index_tmp(i) = i;
    f_tmp(i) = Population(i).f;
    phi_tmp(i) = Population(i).phi;
end


%% For debug
% for i = 1 : n_population - 1
%     index1 = randi(n_population);
%     index2 = randi(n_population);
%     [index_tmp(index1), index_tmp(index2)] = swap(index_tmp(index1),index_tmp(index2));
%     [f_tmp(index1), f_tmp(index2)] = swap(f_tmp(index1),f_tmp(index2));
%     [phi_tmp(index1), phi_tmp(index2)] = swap(phi_tmp(index1),phi_tmp(index2));
% end


%% Stochastic ranking sort
for i = 1 : n_population
    flg = 0;
    for j = 1 : n_population - 1
        u = rand;
        if phi_tmp(j) == phi_tmp(j+1) || u < Pf
            if f_tmp(j) > f_tmp(j+1)
                [index_tmp(j), index_tmp(j+1)]  = RCGAswap(index_tmp(j),index_tmp(j+1));
                [f_tmp(j), f_tmp(j+1)] = RCGAswap(f_tmp(j),f_tmp(j+1));
                [phi_tmp(j), phi_tmp(j+1)] = RCGAswap(phi_tmp(j),phi_tmp(j+1));
                flg = 1;
            end
        else
            if phi_tmp(j) > phi_tmp(j+1)
                [index_tmp(j), index_tmp(j+1)] = RCGAswap(index_tmp(j),index_tmp(j+1));
                [f_tmp(j), f_tmp(j+1)] = RCGAswap(f_tmp(j),f_tmp(j+1));
                [phi_tmp(j), phi_tmp(j+1)] = RCGAswap(phi_tmp(j),phi_tmp(j+1));
                flg = 1;
            end
        end
    end
    if flg == 0
        break;
    end
end


%% Making output variable
Sorted_Population(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);

for i = 1 : n_population
    Sorted_Population(i) = Population(index_tmp(i));
end
