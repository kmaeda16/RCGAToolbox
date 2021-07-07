function [] = RCGAToolbox_Install_ALL
% Installation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

diary('RCGAToolbox_Install_Log.txt');


fprintf('################################################################\n');
fprintf('#                RCGAToolbox Installation Tool                 #\n');
fprintf('################################################################\n');
fprintf('\n');


cd ..;
rcgatoolbox_root = pwd;
cd install;


% RCGAToolbox
answ = input('    Install RCGAToolbox core components? (y/n) ','s');
if answ == 'y'
    RCGAToolbox_Install;
end


% IQM Tools
answ = input('    Install IQM Tools? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/IQMtools V1.2.2.2');
    cd(f);
    RCGAToolbox_Install_IQMTools;
end


% SundialsTB
answ = input('    Install SundialsTB? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB');
    cd(f);
    RCGAToolbox_Install_SundialsTB;
end


% libSBML
answ = input('    Install libSBML? (y/n) ','s');
if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/libSBML-5.18.0-matlab-binaries');
    cd(f);
    RCGAToolbox_Install_libSBML;
end


f = fullfile(rcgatoolbox_root,'install');
cd(f);


fprintf('RCGAToolbox installation, all done.\n');


diary off;
