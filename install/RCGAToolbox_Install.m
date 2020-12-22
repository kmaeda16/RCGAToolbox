% Installation script
% This script must be executed under the directory RCGAToolbox/install/
% just after starting MATLAB.

cd ..;
currentdir = pwd;
cd install;
f = fullfile(currentdir,'source');
addpath(genpath(f));
savepath;

fprintf('Done!\n');
