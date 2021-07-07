function [] = RCGAToolbox_Install_libSBML
% libSBML installation script
% This script must be executed under the directory
% RCGAToolbox/3rdparty/libSBML-5.18.0-matlab-binaries/
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

addpath(pwd);
savepath;

fprintf('libSBML was added to the MATLAB path, and the current path settings were saved.\n');
