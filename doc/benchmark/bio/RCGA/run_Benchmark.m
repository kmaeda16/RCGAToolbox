clear all;

n_repeat = 5;


%% Start Calculation
cwd = pwd;
addpath(cwd);
dirname = 'Results';
mkdir(dirname);
cd(dirname);

for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
%         Benchmark_UNDXMGG(name,i); % For normal calculation
        batch(@Benchmark_UNDXMGG,0,{name,i}); % For batch calculation (Parallel Computing Toolbox required)
%         Benchmark_REXstarJGG(name,i); % For normal calculation
        batch(@Benchmark_REXstarJGG,0,{name,i}); % For batch calculation (Parallel Computing Toolbox required)
    end
end

cd('..');
rmpath(cwd);
