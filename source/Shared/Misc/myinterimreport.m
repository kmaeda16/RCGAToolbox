function interimreport(elapsedTime,i,Param,best)

printTransition(elapsedTime,i,best);
writeTransition(elapsedTime,i,Param,best);

if ~isfield(Param,'interimreportfun')
    return;
end
interimreportfun = Param.interimreportfun;
decodingfun = Param.decodingfun;
n_constraint = Param.n_constraint;

x = decodingfun(best.gene);

if n_constraint == 0
    interimreportfun(i,x,best.f);
else
    interimreportfun(i,x,best.f,best.phi,best.g);
end
