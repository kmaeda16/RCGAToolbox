function [] = RCGAToolbox_Uninstall_ALL
% Uninstallation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

cd ..;
currentdir = pwd;
cd install;

% RCGAToolbox
f = fullfile(currentdir,'source');
rmpath(genpath(f));

% IQM Tools
f = fullfile(currentdir,'3rdparty/IQMtools V1.2.2.2');
rmpath(genpath(f));

% SundialsTB
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/cvodes');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/cvodes/cvm');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/cvodes/function_types');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/idas');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/idas/idm');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/idas/function_types');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/kim');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/function_types');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/nvector');
rmpath(f);
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB/putils');
rmpath(f);

% libSBML
f = fullfile(currentdir,'3rdparty/libSBML-5.18.0-matlab-binaries');
rmpath(genpath(f));

savepath;

fprintf('RCGAToolbox and the third-party software tools were removed from the MATLAB path, and the current path settings were saved.\n');
