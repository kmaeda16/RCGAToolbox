function RCGAcreateMeasurement(PEtabMeasurementFile, IQMmeasurementFile, model)
% RCGAcreateMeasurement creates a IQM Tools measurement file from PEtab 
% measurement file.
% 
% [SYNTAX]
% RCGAcreateMeasurement(PEtabMeasurementFile, IQMmeasurementFile, model)
% 
% [INPUT]
% PEtabMeasurementFile :  Name of PEtab measurement file. This is a TSV 
%                         file with at least 4 fields: observableId,
%                         simulationConditionId, measurement, and time. For 
%                         the format, see Schmiester, PloS Comput Biol, 
%                         2021, DOI: 10.1371/journal.pcbi.1008646
% IQMmeasurementFile   :  Name of IQM Tools measurement file to be created
%                         from the above PEtab measurement file.
% model                :  An IQMmodel object, the name of SBML file, the 
%                         function handle for an ODE function (IQM Tools 
%                         format), or the function handle for a MEXed 
%                         model.


%% Checking input arguments
if nargin ~= 3
    error('Incorrect number of input arguments.');
end


%% Checking if IQM Tools are available.
if exist('isIQMmodel','file') == 0 || exist('IQMmodel','file') == 0
    warning('IQM Tools are not properly installed. Run the script RCGAToolbox/install/RCGAToolbox_Diagnosis for diagnosis.');
end


%% Preparation
% PEtabMeasurementFile
T = tdfread(PEtabMeasurementFile);
if isfield(T,'simulationConditionId')
    warning('Note that simulationConditionId is not yet supported. It was assumed that all measurement data in %s were obtained in a single experiment.',PEtabMeasurementFile);
end

% IQMmeasurementFile
[ ~, funcname, ~ ] = fileparts(IQMmeasurementFile);
IQMmeasurementFile = strcat(funcname,'.csv');

% model
if isa(model,'function_handle')
    statename = model('states');
    n_state = length(statename);
    
elseif isIQMmodel(model)
    
    st_model = struct(model);
    n_state = length(st_model.states);
    statename = cell(1,n_state);
    for i = 1 : n_state
        statename(i) = { st_model.states(i).name };
    end
    
elseif ischar(model)
    
    [~, name, ext] = fileparts(model);
    
    if strcmp(ext,'.xml') || strcmp(ext,'.sbml')
        model = IQMmodel(model);
        st_model = struct(model);
        n_state = length(st_model.states);
        statename = cell(1,n_state);
        for i = 1 : n_state
            statename(i) = { st_model.states(i).name };
        end
        
    elseif strcmp(ext,'.m') || strcmp(ext(2:end),mexext)
        statename = eval( sprintf('%s(''states'')',name) );
        n_state = length(statename);
        
    else
        error('Unexpected file extension for model!');
        
    end
    
else
    error('Unexpected input for model!');
    
end


%% Making Values matrix for the output IQM measurement file

% Number of data points in the input PEtab measurement file
temp = size(T.observableId);
n_datapoint = temp(1);

% Number of unique time points in the input PEtab measurement file
uniqtime = unique(T.time);
n_timepoint = length(uniqtime);

Values = [ uniqtime, nan(n_timepoint,n_state) ]; 

for i = 1 : n_datapoint
    for j = 1 : n_timepoint
        if T.time(i) == Values(j,1)
            for k = 1 : n_state
                if strcmp(strtrim(T.observableId(i,:)),statename(k))
                    Values(j,k+1) = T.measurement(i);
                end
            end
        end
    end
end


%% Opening Output File
out = fopen(IQMmeasurementFile,'w');
if out == -1
    warning('cannot open %s!\n',IQMmeasurementFile);
    return;
end

%% [Name]
fprintf(out,'[Name]\n');
fprintf(out,'%s\n',PEtabMeasurementFile);
fprintf(out,'\n');

%% [Notes]
fprintf(out,'[Notes]\n');
fprintf(out,'This function was created by RCGAToolbox RCGAcreateMeasurement\n');
fprintf(out,'based on %s.\n',PEtabMeasurementFile);
fprintf(out,'Created: %s\n\n',date);


%% [Components]
fprintf(out,'[Components]\n');
fprintf(out,'time,');
for i = 1 : n_state
    if i < n_state
        fprintf(out,'%s,',char(statename(i)));
    else
        fprintf(out,'%s\n',char(statename(i)));
    end
end
fprintf(out,'\n');


%% [Componentnotes]
fprintf(out,'[Componentnotes]\n');
fprintf(out,'\n');


%% [Values]
fprintf(out,'[Values]\n');
[n_row, n_col] = size(Values);
for i = 1 : n_row
    for j = 1 : n_col
        if j < n_col
            fprintf(out,'%e,',Values(i,j));
        else
            fprintf(out,'%e\n',Values(i,j));
        end
    end
end


%% Print Path and File name
fprintf('Measurement data file "%s" was created in "%s".\n',IQMmeasurementFile,pwd);


%% Closing Output File
fclose(out);
