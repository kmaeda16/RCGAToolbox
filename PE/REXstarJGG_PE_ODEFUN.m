function x = REXstarJGG_PE_ODEFUN(odefun,IC,mst,fast_flag,Param,fitnessfun_PE)

if ~isSBmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = SBmeasurement(mst);
    fprintf(' Finished.\n');
end

if isfield(Param,'opts')
    opts = Param.opts;
else
    opts = struct;
end
Param.fitnessfun = @(x) fitnessfun_PE(x,odefun,IC,mst,fast_flag,opts);
interimreportfun = Param.interimreportfun;
Param.interimreportfun = @(elapsedTime,generation,Param,Population,best) interimreportfun(elapsedTime,generation,Param,Population,best,odefun,IC,mst,fast_flag);

Param = checkInputs(Param,'REXstarJGG');

best = RCGA_Main(Param,@JGG);

x = Param.decodingfun(best.gene);
