% IQM Tools installation script
% This script must be executed under the directory
% RCGAToolbox/3rdparty/IQMtools V1.2.2.2/
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

flg = 0;


try
    
    archstr = computer('arch');
    if strcmp(archstr,'maci64')
        movefile('IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25.c',...
            'IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25_temp.c');
        movefile('IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25_Mac.c',...
            'IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25.c');
    end
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
