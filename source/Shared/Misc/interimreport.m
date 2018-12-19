function interimreport(i,chrom,Param)

interimreportfun = Param.interimreportfun;
decodingfun = Param.decodingfun;
n_constraint = Param.n_constraint;

x = decodingfun(chrom.gene);

if n_constraint == 0
    interimreportfun(i,x,chrom.f);
else
    interimreportfun(i,x,chrom.f,chrom.phi,chrom.g);
end
