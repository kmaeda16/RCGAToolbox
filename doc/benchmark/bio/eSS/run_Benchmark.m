
n_repeat = 5;


%% Start Calculation
for Name = {'hiv','threestep'}
    name = char(Name);
    for i = 1 : n_repeat
        dirname = sprintf('result_%s_%d',name,i);
        mkdir(dirname);
        cd(dirname);
        addpath('..');
%         Benchmark_eSS(name,i); % For normal calculation
        batch(@Benchmark_eSS,0,{name,i}); % For batch calculation
        rmpath('..');
        cd('..');
    end
end


%% Make Result files
% After finishing calculation, execute below in the directory eSS.
% It takes a while.

% dirname= 'Results';
% mkdir(dirname);
% cd(dirname);
% 
% n_repeat = 5;
% 
% 
% addpath(genpath('../../function'));
% 
% for Name = {'hiv','threestep'}
%     name = char(Name);
%     for i = 1 : n_repeat
%         fprintf('%s %d ...\n',name,i);
%         infilename = sprintf('../result_%s_%d/ess_report.mat',name,i);
%         outfilename = sprintf('eSS_%s_transition_%d.dat',name,i);
%         switch name
%             case 'hiv'
%                 makeResultFiles(infilename,outfilename,@wrapper_hiv_con_mex);
%             case 'threestep'
%                 makeResultFiles(infilename,outfilename,@wrapper_threestep_con_mex);
%             otherwise
%                 error('Unexpected problem name.');
%         end  
%     end
% end
% 
% rmpath(genpath('../../function'));
% 
% cd('..');
