function Results = RCGA_REXstarJGG(problem,varargin)
% RCGA_REXstarJGG is a function of REXstar/JGG
% 
% [SYNTAX]
% Results = RCGA_REXstarJGG(problem)
% Results = RCGA_REXstarJGG(problem,opts)
% 
% [INPUT]
% problem :  Problem structure.
% opts    :  RCGA options. See XXXXXXXXXXX for options.
% 
% [OUTPUT]
% Results :  Structure with results


%% Handling inputs
switch nargin
    case 1
        opts = struct;
    case 2
        opts = varargin{1};
    otherwise
        error('Incorrect number of input arguments');
end


%% Checking inputs
[problem, opts] = RCGAcheckInputs(problem,opts,mfilename);


%% Printing welcome messages
RCGAprintWelcomeMessage(problem,opts,mfilename);


%% Executing REXstar/JGG
Results = RCGA_Main(problem,opts,@RCGA_JGG);
