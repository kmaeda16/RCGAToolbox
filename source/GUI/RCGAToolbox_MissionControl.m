function RCGAToolbox_MissionControl(matlab_version)
% RCGAToolbox_MissionControl starts GUI for general optimization problem.
% 
% [SYNTAX]
% RCGAToolbox_MissionControl
% RCGAToolbox_MissionControl(matlab_version)
% 
% [INPUT]
% matlab_version : String that shows MATLAB version, e.g. 'R2016a' or 
%                  'R2020b'. Based on matlab_version, 
%                  RCGAToolbox_MissionControl selects the latest available 
%                  version of GUI. If matlab_version is not specified, 
%                  RCGAToolbox_MissionControl automatically checks the 
%                  MATLAB version on your system.


flg = 0;


if ~exist('matlab_version','var')
    matlab_version = [ 'R' version('-release') ];
end


if lower(matlab_version(end)) == 'a' && lower(matlab_version(1)) == 'r' && length(matlab_version) == 6
    
    matlab_version_num = str2double(matlab_version(2:end-1));
    
elseif lower(matlab_version(end)) == 'b' && lower(matlab_version(1)) == 'r' && length(matlab_version) == 6
    
    matlab_version_num = str2double(matlab_version(2:end-1)) + 0.5;
    
else
    
    matlab_version_num = nan;
    flg = 1;
    
end


if isnan(matlab_version_num)
    flg = 1;
end


if flg == 1
    error('Invalid MATLAB Version. matlab_version should be something like ''R2016a'' or ''R2020b''.');
end


% Newer than or equal to R2021a
if matlab_version_num >= 2021
    
    RCGAToolbox_MissionControl_R2021a;
    
% Older than R2021a
else
    
    RCGAToolbox_MissionControl_R2016a;
    
end
