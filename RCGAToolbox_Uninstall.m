% Uninstallation script

currentdir = pwd;
f = fullfile(currentdir,'source');
rmpath(genpath(f));
savepath;
