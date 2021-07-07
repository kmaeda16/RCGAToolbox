function [] = RCGAToolbox_Uninstall_ALL
% Uninstallation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

diary('RCGAToolbox_Uninstall_Log.txt');


fprintf('################################################################\n');
fprintf('#               RCGAToolbox Uninstallation Tool                #\n');
fprintf('################################################################\n');
fprintf('\n');


cd ..;
rcgatoolbox_root = pwd;
cd install;


% RCGAToolbox
answ = input('    Uninstall RCGAToolbox core components? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'source');
    rmpath(genpath(f));
    savepath;
    fprintf('RCGAToolbox core components were removed from the MATLAB path.\n');
end


% IQM Tools
answ = input('    Uninstall IQM Tools? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/IQMtools V1.2.2.2');
    rmpath(genpath(f));
    savepath;
    fprintf('IQM Tools were removed from the MATLAB path.\n');
end


% SundialsTB
answ = input('    Uninstall SundialsTB? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/cvodes');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/cvodes/cvm');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/cvodes/function_types');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/idas');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/idas/idm');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/idas/function_types');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/kim');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/kinsol/function_types');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/nvector');
    rmpath(f);
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB/putils');
    rmpath(f);
    savepath;
    fprintf('SundialsTB was removed from the MATLAB path.\n');
end


% libSBML
answ = input('    Uninstall libSBML? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/libSBML-5.18.0-matlab-binaries');
    rmpath(genpath(f));
    savepath;
    fprintf('libSBML was removed from the MATLAB path.\n');
end


fprintf('RCGAToolbox uninstallation, all done.\n');


diary off;
