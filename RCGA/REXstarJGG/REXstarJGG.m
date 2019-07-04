function Results = REXstarJGG(problem,varargin)
% REXstarJGG is a function of REXstar/JGG
% 
% [SYNTAX]
% Results = REXstarJGG(problem)
% Results = REXstarJGG(problem,opts)
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
[problem, opts] = checkInputs(problem,opts,mfilename);


%% Printint welcome messages
printWelcomeMessage(problem,opts,mfilename);


%% Executing REXstar/JGG
Results = RCGA_Main(problem,opts,@JGG);
