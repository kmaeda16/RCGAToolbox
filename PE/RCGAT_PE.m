function RCGAT_PE(model,expdata,RCGAParam)

Param.n_population = 200;
Param.n_children = 200;
Param.n_generation = 1e+5;
Param.n_parent = Param.n_gene + 1;
Param.t_rextar = 6.0;
Param.output_intvl = 1e+8;
Param.selection_type = 0;
Param.t_limit = 24 * 60 * 60;
Param.par = 0;
Param.out_population = 'None'; % 'Population.dat';
Param.out_solution = 'None'; % 'Solution.dat';
Param.fitnessfun = @SSR;
Param.decodingfun = mydecodingfun;
Param.model = model;
Param.expdata = expdata;


[ best, Population ] = REXstarJGG(Param);

% input = model, expdata

fprintf('%s\n',char(Problem_Name));
Param = getParam(Param,char(Problem_Name));



SBsimulate