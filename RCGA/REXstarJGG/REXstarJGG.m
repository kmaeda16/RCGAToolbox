function [ best, Population ] = REXstarJGG(problem,varargin)

switch nargin
    case 1
        opts = struct;
    case 2
        opts = varargin{1};
    otherwise
        error('Incorrect number of input arguments');
end

[problem, opts] = checkInputs(problem,opts,mfilename);

[ best, Population ] = RCGA_Main(problem,opts,@JGG);
