% Uninstallation script for RCGAToolbox only (no third-party software tools uninstalled)
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

cd ..;
currentdir = pwd;
cd install;
f = fullfile(currentdir,'source');
rmpath(genpath(f));
savepath;

fprintf('RCGAToolbox was removed from the MATLAB path, and the current path settings were saved.\n');
