% Installation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

RCGAToolbox_Install;

cd('../3rdparty/IQMtools V1.2.2.2'); % MEX command and a C compiler compatible with MATLAB required.
RCGAToolbox_Install_IQMTools;

cd('../sundials-2.6.2/sundialsTB');
RCGAToolbox_Install_SundialsTB; % MEX command and a C compiler compatible with MATLAB required.

cd('../libSBML-5.18.0-matlab-binaries');
RCGAToolbox_Install_libSBML;

cd('../../install');
