function [] = RCGAToolbox_Install_ALL
% Installation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

cd ..;
currentdir = pwd;
cd install;

% RCGAToolbox
RCGAToolbox_Install;

% IQM Tools
f = fullfile(currentdir,'3rdparty/IQMtools V1.2.2.2');
cd(f);
RCGAToolbox_Install_IQMTools; % MEX command and a C compiler compatible with MATLAB required.

% SundialsTB
f = fullfile(currentdir,'3rdparty/sundials-2.6.2/sundialsTB');
cd(f);
RCGAToolbox_Install_SundialsTB; % MEX command and a C compiler compatible with MATLAB required.

% libSBML
f = fullfile(currentdir,'3rdparty/libSBML-5.18.0-matlab-binaries');
cd(f);
RCGAToolbox_Install_libSBML;

cd('../../install');
