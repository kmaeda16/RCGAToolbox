% SundialsTB installation script
% This script must be executed under the directory
% RCGAToolbox/3rdparty/sundials-2.6.2/sundialsTB/
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

flg = 0;

try
    install_STB_mod;
    startup_STB_mod;
    savepath;
catch ME
    warning(ME.message);
    flg = 1;
end


if flg == 0
    fprintf('SundialsTB was successfully installed for RCGAToolbox, and the current path settings were saved.\n');
else
    error('Installation of SundialsTB failed.');
end
