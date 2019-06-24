function [ Results, optimizedmodel ] = REXstarJGG_c(model,decodingfun,mst,varargin)
% REXstarJGG_c(model,decodingfun,mst);
% REXstarJGG_c(model,decodingfun,mst,opts); 4
% REXstarJGG_c(model,decodingfun,mst,simopts,opts); 5 
% REXstarJGG_c(model,decodingfun,mst,fast_flag,simopts,opts); 6
% REXstarJGG_c(model,decodingfun,mst,n_gene,n_constraint,fitnessfun,fast_flag,simopts,opts);

switch nargin
    case 3
        fast_flag = 1;
        fitnessfun = @SSR_c;
        n_constraint = 0;
        simopts = struct;
        opts = struct;
    case 4
        fast_flag = 1;
        fitnessfun = @SSR_c;
        n_constraint = 0;
        simopts = struct;
        opts = varargin{1};
    case 5
        fast_flag = 1;
        fitnessfun = @SSR_c;
        n_constraint = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 6
        fast_flag = varargin{1};
        fitnessfun = @SSR_c;
        n_constraint = 0;
        simopts = varargin{2};
        opts = varargin{3};
    case 9
        n_gene = varargin{1};
        n_constraint = varargin{2};
        fitnessfun = varargin{3};
        fast_flag = varargin{4};
        simopts = varargin{5};
        opts = varargin{6};
    otherwise
        error('Incorrect number of input arguments');
end


if isempty(n_constraint)
    n_constraint = 0;
end
if isempty(fitnessfun)
    fitnessfun = @SSR_c;
end
if isempty(fast_flag)
    fast_flag = 1;
end
if isempty(simopts)
    simopts = struct;
end
if isempty(opts)
    opts = struct;
end


% If model is an SBmodel object and fast_flg is one
fprintf('Compiling %s.c with %s.h...\n',model,model);
mexcompileIQM(model);
mex_name = model;


if isempty(which('IQMmodel'))
    error('IQMTools is not properly installed.');
end

if ~isIQMmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = IQMmeasurement(mst);
    fprintf(' Finished.\n');
end
st_mst = struct(mst{1});
n_var = length(st_mst.data);



problem.n_gene = n_gene;
problem.n_constraint = n_constraint;
problem.fitnessfun = @(x) fitnessfun(x,mst,mex_name,n_var,simopts);
problem.decodingfun = decodingfun;

if ~isfield(opts,'interimreportfun')
    opts.interimreportfun = @interimreportfun_c;
end
interimreportfun = opts.interimreportfun;
opts.interimreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
    interimreportfun(...
    elapsedTime,generation,problem,opts,Population,best,...
    mst,mex_name,n_var,simopts,fast_flag);

Results = REXstarJGG(problem,opts);

