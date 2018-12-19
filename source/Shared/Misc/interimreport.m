function interimreport(chrom,Param)

interimreportfun = Param.interimreportfun;
decodingfun = Param.decodingfun;
n_constraint = Param.n_constraint;

x = decodingfun(chrom.gene);

if n_constraint == 0
    interimreportfun(x,chrom.f);
else
    interimreportfun(x,chrom.f,chrom.phi,chrom.g);
end
