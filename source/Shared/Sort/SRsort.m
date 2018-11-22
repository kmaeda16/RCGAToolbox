function Population = SRsort(Population, Pf)

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
	
for i = 1 : n_population - 1
    index1 = randi(n_population);
    index2 = randi(n_population);
    [index_tmp(index1), index_tmp(index2)] = swap(index_tmp(index1),index_tmp(index2));
    [f_tmp(index1), f_tmp(index2)] = swap(f_tmp(index1),f_tmp(index2));
    [phi_tmp(index1), phi_tmp(index2)] = swap(phi_tmp(index1),phi_tmp(index2));
end
	
for i = 1 : n_population
    flg = 0;
    for j = 1 : n_population - 1
        u = rand;
        if phi_tmp(j) == phi_tmp(j+1) || u < Pf
            if f_tmp(j) > f_tmp(j+1)
               [index_tmp(j), index_tmp(j+1)]  = swap(index_tmp(j),index_tmp(j+1));
					[f_tmp(j), f_tmp(j+1)] = swap(f_tmp(j),f_tmp(j+1));
					[phi_tmp(j),phi_tmp(j+1)] = swap(phi_tmp(j),phi_tmp(j+1));
					flg = 1;
            end
        else
            if phi_tmp(j) >  phi_tmp(j+1)
                [index_tmp(j), index_tmp(j+1)] = swap(index_tmp(j),index_tmp(j+1));
                [f_tmp(j), f_tmp(j+1)] = swap(f_tmp(j),f_tmp(j+1));
                [phi_tmp(j),phi_tmp(j+1)] = swap(phi_tmp(j),phi_tmp(j+1));
                flg = 1;
            end
        end
    end
    if flg == 0
        break;
    end
end

pop_tmp(1,1:n_population) = struct('gene',zeros(1,n_gene),'g',zeros(1,n_constraint),'f',0,'phi',0);

for i = 1 : n_population
    pop_tmp(i) = Population(i);
end
for i = 1 : n_population
    Population(i) = pop_tmp(index_tmp(i));
end
