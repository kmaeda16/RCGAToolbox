% Uninstallation script
% This script must be executed under the directory RCGAToolbox/install/
% just after starting MATLAB.

cd ..;
currentdir = pwd;
cd install;
f = fullfile(currentdir,'source');
rmpath(genpath(f));
savepath;

fprintf('Done!\n');
