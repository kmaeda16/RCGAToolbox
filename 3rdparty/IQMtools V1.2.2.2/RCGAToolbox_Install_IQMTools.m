% IQM Tools installation script
% This script must be executed under the directory
% RCGAToolbox/3rdparty/IQMtools V1.2.2.2/
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

flg = 0;

try
    installIQMtoolsInitial;
    savepath;
catch ME
    warning(ME.message);
    flg = 1;
end


if flg == 0
    fprintf('IQM Tools were successfully installed, and the current path settings were saved.\n');
else
    error('Installation of IQM Tools failed.');
end
