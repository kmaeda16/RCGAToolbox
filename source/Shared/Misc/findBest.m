function index = findBest(Population)

n_population = length(Population);
f = Inf;
phi = Inf;
index = 0;

for i = 1 : n_population
    if Population(i).phi < phi ...
            || ( Population(i).phi == phi && Population(i).f < f )
        f = Population(i).f;
        phi = Population(i).phi;
        index = i;
    end
end
