function mexfilename = RCGAmakeMEXmodel(model, mexfilename, doNOTcompileFlag)
% RCGAmakeMEXmodel makes a MEX file based on model. Note that a C source
% code (*.c) and header file (*.h) are always generated.
% 
% [SYNTAX]
% mexfilename = RCGAmakeMEXmodel(model)
% mexfilename = RCGAmakeMEXmodel(model, mexfilename)
% mexfilename = RCGAmakeMEXmodel(model, mexfilename, doNOTcompileFlag)
% 
% [INPUT]
% model            :  An IQM object or the name of a SBML file (*.sbml or 
%                     *.xml).
% mexfilename      :  Name of the generated MEX file.
% doNOTcompileFlag :  Flag specifies if a MEX file is compiled (generated).
%                     - doNOTcompileFlag = 0 : MEX file is compiled
%                        (default).
%                     - doNOTcompileFlag = 1 : MEX file is NOT compiled.
% 
% [OUTPUT]
% mexfilename      :  Name of the generated MEX file.


%% Checking if IQM Tools are available.
if exist('isIQMmodel','file') == 0 || ...
        exist('IQMmakeMEXmodel','file') == 0 || ...
        exist('mexcompileIQM','file') == 0 || ...
        exist('IQMmodel','file') == 0
    warning('IQM Tools are not properly installed. Run the script RCGAToolbox/install/RCGAToolbox_Diagnosis for diagnosis.');
end


%%
if ~exist('doNOTcompileFlag','var')
    doNOTcompileFlag = 0;
end

if isIQMmodel(model)
    sbm = model;
    if ~exist('mexfilename','var')
        st_sbm = struct(sbm);
        mexfilename = strcat(st_sbm.name,'_mex');
        mexfilename = regexprep(mexfilename,'\W','');
    end
elseif ischar(model)
    sbm = IQMmodel(model);
    if ~exist('mexfilename','var')
        [ ~, filename, ~] = fileparts(model);
        mexfilename = strcat(filename,'_mex');
    end
else
    error('model must be an IQMmodel or a file name!');
end


IQMmakeMEXmodel(sbm,mexfilename,1);

if doNOTcompileFlag == 0
    mexcompileIQM(mexfilename);
elseif doNOTcompileFlag ~= 1
    error('Unexpected doNOTcompileFlag!');
end
