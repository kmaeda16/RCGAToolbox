function [] = RCGAToolbox_Install
% Installation script for RCGAToolbox only (no third-party software tools installed)
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

cd ..;
currentdir = pwd;
cd install;
f = fullfile(currentdir,'source');
addpath(genpath(f));
savepath;

fprintf('RCGAToolbox was added to the MATLAB path, and the current path settings were saved.\n');
