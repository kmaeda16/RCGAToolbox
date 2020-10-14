function Results = RCGA_CustomRCGA(problem, GenerationAlteration, varargin)
% RCGA_CustomRCGA is a function that helps users to implement their own
% real-coded genetic algorithms (RCGAs). All users have to do is to define
% a generation alteration function. Use
% RCGAToolbox/doc/demo/CustomRCGA/ExampleGenerationAlteration.m as a
% template.
% 
% 
% [SYNTAX]
% Results = RCGA_CustomRCGA(problem,GenerationAlteration)
% Results = RCGA_CustomRCGA(problem,GenerationAlteration,opts)
% 
% [INPUT]
% problem              :  Problem structure.
% GenerationAlteration :  Function handle for custom generation alteration 
%                         function.
% opts                 :  RCGA options. See XXXXXXXXXXX for options.
% 
% [OUTPUT]
% Results :  Structure with results
% 


%% Handling inputs
switch nargin
    case 2
        opts = struct;
    case 3
        opts = varargin{1};
    otherwise
        error('Incorrect number of input arguments');
end


%% Checking inputs
[problem, opts] = RCGAcheckInputs(problem,opts,mfilename);


%% Printing welcome messages
RCGAprintWelcomeMessage(problem,opts,mfilename);


%% Executing CustomRCGA
Results = RCGA_Main(problem,opts,GenerationAlteration);
