
n_repeat = 5;
dirname = 'Results';

mkdir(dirname);
cd(dirname);
addpath('..');

for Name = {'hiv','threestep'}
    name = char(Name);
    
    for n_par = [ 1, 2, 4, 8, 16, 32 ]
        
        if 1 < n_par
            p = gcp('nocreate');
            if isempty(p)
                parpool(n_par);
            elseif ~isempty(p)
                if p.NumWorkers ~= n_par
                    delete(p);
                    parpool(n_par);
                end
            end
        end
        
        for i = 1 : n_repeat
            Benchmark_UNDXMGG(name,i,n_par); % For normal calculation
            Benchmark_REXstarJGG(name,i,n_par); % For normal calculation
        end
        
    end
    
end

rmpath('..');
cd('..');
