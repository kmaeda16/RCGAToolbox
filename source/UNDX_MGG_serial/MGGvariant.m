function population = MGGvariant(population, n_children, SearchRegion)
% Minimal Generation Gap Selection for UNDX
% See Hiroaki Kitano, "Genetic Algorithms 4", Sangyo-tosho, p234, 2000 

maxitr = 10;

n_population = length(population);

ip = zeros(1,3);
for i = 1:maxitr
    ip(1)=ceil(rand*n_population);
    ip(2)=ceil(rand*n_population);
    if ip(1) ~= ip(2)
        break;
    end
end

for i = 1:n_children
    for j = 1:maxitr
        ip(3) = ceil(rand*n_population);
        if ip(3) ~= ip(1) && ip(3) ~= ip(2)
            break;
        end
    end
    children(i) = getNewChild(population(ip(1)),population(ip(2)),population(ip(3)),SearchRegion);
end

children(n_children+1) = population(ip(1));
children(n_children+2) = population(ip(2));
children = mysort(children);
population(ip(1)) = children(1);
population(ip(2)) = children(ceil(rand*(n_children+1))+1);
population = mysort(population);
