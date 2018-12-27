function x = REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene,varargin)
% REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene); 5
% REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene,opts); 6
% REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene,simopts,opts); 7
% REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene,fast_flag,simopts,opts); 8
% REXstarJGG_odefun(odefun,icfun,decodingfun,mst,n_gene,n_constraint,fitnessfun,fast_flag,simopts,opts); 10

switch nargin
    case 5
        fast_flag = 1;
        fitnessfun = @SSR_odefun;
        n_constraint = 0;
        simopts = struct;
        opts = struct;
    case 6
        fast_flag = 1;
        fitnessfun = @SSR_odefun;
        n_constraint = 0;
        simopts = struct;
        opts = varargin{1};
    case 7
        fast_flag = 1;
        fitnessfun = @SSR_odefun;
        n_constraint = 0;
        simopts = varargin{1};
        opts = varargin{2};
    case 8
        fast_flag = varargin{1};
        fitnessfun = @SSR_odefun;
        n_constraint = 0;
        simopts = varargin{2};
        opts = varargin{3};
    case 10 
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
    fitnessfun = @SSR_odefun;
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

% fast_flag = 0;
if fast_flag == 1 && isempty(strfind(which('cvm'),'sundialsTB'))
    warning('sundialsTB is not properly installed. fast_flag set to 0.');
end

if isempty(which('IQMmodel'))
    error('IQMTools is not properly installed.');
end

if ~isIQMmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = IQMmeasurement(mst);
    fprintf(' Finished.\n');
end


problem.n_gene = n_gene;
problem.n_constraint = n_constraint;
problem.fitnessfun = @(x) fitnessfun(x,odefun,icfun,mst,fast_flag,simopts);
problem.decodingfun = decodingfun;

if ~isfield(opts,'interimreportfun')
    opts.interimreportfun = @interimreportfun_odefun;
end
interimreportfun = opts.interimreportfun;
opts.interimreportfun = @(elapsedTime,generation,problem,opts,Population,best) ...
    interimreportfun(...
    elapsedTime,generation,problem,opts,Population,best,...
    odefun,icfun,mst,simopts,fast_flag);

best = REXstarJGG(problem,opts);

x = decodingfun(best.gene);
