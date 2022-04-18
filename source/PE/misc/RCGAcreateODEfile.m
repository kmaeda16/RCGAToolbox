function odefilename = RCGAcreateODEfile(model, odefilename)
% RCGAcreateODEfile reads model (an IQMmodel object or an SBML file) and
% creates an ODE file (IQM Tools format).
% 
% [SYNTAX]
% odefilename = RCGAcreateODEfile(model)
% odefilename = RCGAcreateODEfile(model, odefilename)
% 
% [INPUT]
% model       :  An IQM object or the name of a SBML file (*.sbml or 
%                *.xml).
% odefilename :  Name of the created ODE file (IQM Tools format).
% 
% [OUTPUT]
% odefilename :  Name of the created ODE file (IQM Tools format).


%% Checking if IQM Tools are available.
if exist('isIQMmodel','file') == 0 || ...
        exist('IQMcreateODEfile','file') == 0 || ...
        exist('IQMmodel','file') == 0
    warning('IQM Tools are not properly installed. Run the script RCGAToolbox/install/RCGAToolbox_Diagnosis for diagnosis.');
end


%%
if isIQMmodel(model)
    sbm = model;
    if ~exist('odefilename','var')
        st_sbm = struct(sbm);
        odefilename = strcat(st_sbm.name,'_odefun');
        odefilename = regexprep(odefilename,'\W','');
    end
elseif ischar(model)
    sbm = IQMmodel(model);
    if ~exist('odefilename','var')
        [ filepath, filename, ~] = fileparts(model);
        odefilename = strcat(fullfile(filepath,filename),'_odefun');
    end
else
    error('model must be an IQMmodel or a file name!');
end

IQMcreateODEfile(sbm,odefilename);
