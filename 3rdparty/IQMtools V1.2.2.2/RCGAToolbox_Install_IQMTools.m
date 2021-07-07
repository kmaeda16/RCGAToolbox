function [] = RCGAToolbox_Install_IQMTools
% IQM Tools installation script
% This script must be executed under the directory
% RCGAToolbox/3rdparty/IQMtools V1.2.2.2/
% Run this script just after starting MATLAB. Otherwise, unnecessary paths 
% might be saved.

mexcompiler = 'mex -v';

mex_ok = check_mex(mexcompiler);

if ~mex_ok
  return
end


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
