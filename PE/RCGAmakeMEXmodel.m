function mexfilename = RCGAmakeMEXmodel(model, mexfilename, doNOTcompileFlag)
% RCGAmakeMEXmodel makes a MEX file based on model. Note that C source
% code is always generated.
% 
% [SYNTAX]
% mexfilename = RCGAmakeMEXmodel(model)
% mexfilename = RCGAmakeMEXmodel(model, mexfilename)
% mexfilename = RCGAmakeMEXmodel(model, mexfilename, doNOTcompileFlag)
% 
% [INPUT]
% model            : IQMmodel or SBML file name (*.sbml or *.xml)
% mexfilename      : Name of MEX file
% doNOTcompileFlag : Flag specifies if a MEX file is compiled (generated)
%                    * doNOTcompileFlag = 0 : MEX file is compiled (default)
%                    * doNOTcompileFlag = 1 : MEX file is NOT compiled
% 
% [OUTPUT]
% mexfilename      : Name of the generated MEX file


if ~exist('doNOTcompileFlag','var')
    doNOTcompileFlag = 0;
end

if isIQMmodel(model)
    sbm = model;
    if ~exist('odefilename','var')
        st_sbm = struct(sbm);
        mexfilename = strcat(st_sbm.name,'_mex');
        mexfilename = regexprep(mexfilename,'\W','');
    end
elseif ischar(model)
    sbm = IQMmodel(model);
    if ~exist('odefilename','var')
        [ ~, filename, ~] = fileparts(model);
        mexfilename = strcat(filename,'_mex');
    end
else
    error('model must be an IQMmodel or file name!');
end


IQMmakeMEXmodel(sbm,mexfilename,1);

if doNOTcompileFlag == 1
    mexcompileIQM(mexfilename);
elseif doNOTcompileFlag ~= 0
    error('Unexpected doNOTcompileFlag!');
end
