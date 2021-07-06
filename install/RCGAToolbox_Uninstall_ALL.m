% Uninstallation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

cd ..;
currentdir = pwd;
cd install;
f = fullfile(currentdir,'source');
rmpath(genpath(f));
f = fullfile(currentdir,'3rdparty');
rmpath(genpath(f));
savepath;

fprintf('RCGAToolbox and the third-party software tools were removed from the MATLAB path, and the current path settings were saved.\n');
