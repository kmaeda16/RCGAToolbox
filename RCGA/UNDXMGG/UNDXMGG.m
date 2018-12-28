function Results = UNDXMGG(problem,varargin)

switch nargin
    case 1
        opts = struct;
    case 2
        opts = varargin{1};
    otherwise
        error('Incorrect number of input arguments');
end

[problem, opts] = checkInputs(problem,opts,mfilename);

Results = RCGA_Main(problem,opts,@MGG);
