function [] = RCGAToolbox_Install
% Installation script for RCGAToolbox, IQM Tools, SundialsTB, and libSBML.
% This script must be executed under the directory RCGAToolbox/install/.
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

diary('RCGAToolbox_Install_Log.txt');


fprintf('################################################################\n');
fprintf('#                RCGAToolbox Installation Tool                 #\n');
fprintf('################################################################\n');
fprintf('Executed: %s\n',char(datetime('now')));
fprintf('\n');


cd ..;
rcgatoolbox_root = pwd;
cd install;


%% RCGAToolbox
answ = input('    Install RCGAToolbox core components? (y/n) ','s');

if answ == 'y'
    f = fullfile(rcgatoolbox_root,'source');
    addpath(genpath(f));
    fprintf('RCGAToolbox core components were added to the MATLAB path.\n\n');
end


%% IQM Tools
answ = input('    Install IQM Tools? (y/n) ','s');

if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/IQMtools V1.2.2.2');
    cd(f);
    flg = 0;
    try
        mexcompiler = 'mex -v';
        mex_ok = check_mex(mexcompiler);
        if mex_ok
            archstr = computer('arch');
            if strcmp(archstr,'maci64')
                movefile('IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25.c',...
                    'IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25_temp.c');
                copyfile(fullfile(rcgatoolbox_root,'install/3rdparty/CVODEmex25_Mac.c'),...
                    'IQMpro/tools/01-MEXmodels/CVODEMEX/src/CVODEmex25.c');
            end
            installIQMtoolsInitial;
        else
            flg = 1;
        end
    catch ME
        warning(ME.message);
        flg = 1;
    end
    if flg == 0
        fprintf('IQM Tools were successfully installed.\n\n');
    else
        warning('Installation of IQM Tools failed.');
    end
end


%% SundialsTB
answ = input('    Install SundialsTB? (y/n) ','s');

if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/sundials-2.6.2/sundialsTB');
    cd(f);
    % flg = 0;
    try
        f = fullfile(rcgatoolbox_root,'install/3rdparty');
        addpath(f);
        flg = install_STB_mod;
        rmpath(f);
    catch ME
        warning(ME.message);
        flg = 1;
    end
    if flg == 0
        copyfile('startup_STB.in','startup_STB.m');
        startup_STB(pwd);
        fprintf('SundialsTB was successfully installed.\n\n');
    else
        warning('Installation of SundialsTB failed.');
    end
end


%% libSBML
answ = input('    Install libSBML? (y/n) ','s');

if answ == 'y'
    f = fullfile(rcgatoolbox_root,'3rdparty/libSBML-5.18.0-matlab-binaries');
    addpath(genpath(f));
    fprintf('libSBML was added to the MATLAB path.\n\n');
    archstr = computer('arch');
    if strcmp(archstr,'maci64')
        warning('When you first use RCGAToolbox SBML-related functions, Mac will block "OutputSBML.mexmaci64" and "TranslateSBML.mexmaci64". Don''t move them to Trash. Go to "Apple menu > System Preferences > Security & Privacy > General", and allow "OutputSBML.mexmaci64" and "TranslateSBML.mexmaci64". ');
        fprintf('\n');
    end
end


%% Finalize
f = fullfile(rcgatoolbox_root,'install');
cd(f);

savepath;
fprintf('The current path settings were saved.\n\n');


fprintf('RCGAToolbox installation, all done.\n\n');


diary off;


%%
%---------------------------------------------------------------------------------
% Check if mex works and if the user accepts the current mexopts
%---------------------------------------------------------------------------------
function mex_ok = check_mex(mexcompiler)

% Create a dummy file
fid = fopen('foo.c', 'w');
fprintf(fid,'#include "mex.h"\n');
fprintf(fid,'void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])\n');
fprintf(fid,'{return;}\n');
fclose(fid);

% Run mexcompiler on foo.c
mex_cmd = sprintf('%s foo.c', mexcompiler);
eval(mex_cmd);

% Remove dummy source file and resulting mex file
delete('foo.c')
delete(sprintf('foo.%s', mexext))

fprintf('\n\nMEX files will be compiled and built using the above options\n');
answ = input('    Proceed? (y/n) ','s');
if answ == 'y'
  mex_ok = true;
else
  fprintf('\n\nOK. All done.\n');
  mex_ok = false;
end

return
