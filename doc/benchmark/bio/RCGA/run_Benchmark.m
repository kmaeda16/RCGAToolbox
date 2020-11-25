
n_repeat = 5;
dirname = 'Results';

mkdir(dirname);
cd(dirname);
addpath('..');

for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
%         Benchmark_UNDXMGG(name,i); % For normal calculation
        batch(@Benchmark_UNDXMGG,0,{name,i}); % For batch calculation
%         Benchmark_REXstarJGG(name,i); % For normal calculation
        batch(@Benchmark_REXstarJGG,0,{name,i}); % For batch calculation
    end
end

rmpath('..');
cd('..');
