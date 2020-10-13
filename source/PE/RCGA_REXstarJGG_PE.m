function Results = RCGA_REXstarJGG_PE(model,decodingfun,mst,varargin)
% RCGA_REXstarJGG_PE estimates parameters included in model by fitting it
% to experimental data mst.
% 
% [SYNTAX]
% Results = RCGA_REXstarJGG_PE(model,decodingfun,mst)
% Results = RCGA_REXstarJGG_PE(model,decodingfun,mst,opts)
% Results = RCGA_REXstarJGG_PE(model,decodingfun,mst,simopts,opts)
% Results = RCGA_REXstarJGG_PE(model,decodingfun,mst,fast_flag,simopts,opts)
% Results = RCGA_REXstarJGG_PE(model,decodingfun,mst,n_constraint, ...
%                         fitnessfun,fast_flag,simopts,opts)
% 
% [INPUT]
% model        :  IQMmodel or file name (*.sbml, *.xml, *.m, *.mex, or *.c)
% decodingfun  :  Function handle for decoding function
% mst          :  Experimental data (IQMmeasurement) or filename
% n_constraint :  Number of constraints.
% fitnessfun   :  Function handle for fitness function
% fast_flag    :  Solver flag
%                 * fast_flag = 0 : ODEXX by MATLAB
%                 * fast_flag = 1 : CVODE by SundialsTB
%                 * fast_flag = 2 : CVODE by IQM Tools
% simopts      :  Structure with integrator options. Fields depend on
%                 Simulation_*. See 'help Simulation_*'.
% opts         :  Structure with RCGA options
% 
% [OUTPUT]
% Results      :  Structure with results


%% Handling input arguments
switch nargin
    case 3
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = struct;
        opts = struct;
    case 4
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = struct;
        opts = varargin{1};
    case 5
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 6
        n_constraint = 0;
        fitnessfun = @RCGAssr;
        fast_flag = varargin{1};
        simopts = varargin{2};
        opts = varargin{3};
    case 8
        n_constraint = varargin{1};
        fitnessfun = varargin{2};
        fast_flag = varargin{3};
        simopts = varargin{4};
        opts = varargin{5};
    otherwise
        error('Incorrect number of input arguments');
end

if isempty(n_constraint)
    n_constraint = 0;
end
if isempty(fitnessfun)
    fitnessfun = @RCGAssr;
end
if isempty(fast_flag)
    fast_flag = 0;
end
if isempty(simopts)
    simopts = struct;
end
if isempty(opts)
    opts = struct;
end


%% Run parameter estimation
Results = RCGA_PE(model,decodingfun,mst,n_constraint,fitnessfun,fast_flag,simopts,opts,@RCGA_REXstarJGG);
