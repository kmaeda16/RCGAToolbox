function RCGAToolbox_MissionControl_PE(matlab_version)
% RCGAToolbox_MissionControl_PE starts GUI for parameter estimation.
% 
% [SYNTAX]
% RCGAToolbox_MissionControl_PE()
% RCGAToolbox_MissionControl_PE(matlab_version)
% 
% [INPUT]
% matlab_version : String that shows MATLAB version, e.g. 'R2016a' or 
%                  'R2020b'. Based on matlab_version, 
%                  RCGAToolbox_MissionControl_PE selects the latest 
%                  available version of GUI. If matlab_version is not 
%                  specified, RCGAToolbox_MissionControl_PE automatically 
%                  checks the MATLAB version on your system.


if ~exist('matlab_version','var')
    matlab_version = [ 'R' version('-release') ];
end


if lower(matlab_version(end)) == 'a'
    
    matlab_version_num = str2double(matlab_version(2:end-1));
    
elseif lower(matlab_version(end)) == 'b'
    
    matlab_version_num = str2double(matlab_version(2:end-1)) + 0.5;
    
else
    
    error('Invalid MATLAB Version. matlab_version should be something like ''R2016a'' or ''R2020b''.');
    
end


% Newer than or equal to R2021a
if matlab_version_num >= 2021
    
    % RCGAToolbox_MissionControl_PE_R2021a;
    RCGAToolbox_MissionControl_PE_R2016a;
    
% Older than R2021a
else
    
    RCGAToolbox_MissionControl_PE_R2016a;
    
end
