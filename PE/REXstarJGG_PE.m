function [ x, optimizedmodel ] = REXstarJGG_PE(model,mst,fast_flg,Param,fitnessfun_PE)

if isSBmodel(model)
    if fast_flg == 1
        temp = SBstruct(model);
        mex_name = strcat(temp.name,'_mex');
        SBPDmakeMEXmodel(model,mex_name);
    else
        mex_name = [];
    end
else
    % If model is the name of a MEX file
    if exist(model,'file') == 3
        mex_name = model;
        model = [];
    else
        % If model is the name of a SBML file
        if exist(model,'file') == 2
            fprintf('Reading %s ...',model);
            model = SBmodel(model);
            fprintf(' Finished.\n');
        end
        % If model is an SBmodel object and fast_flg is one
        if fast_flg == 1
            temp = SBstruct(model);
            mex_name = strcat(temp.name,'_mex');
            SBPDmakeMEXmodel(model,mex_name);
        else
            mex_name = [];
        end
    end
end

Param.fitnessfun = @(x) fitnessfun_PE(x,model,mst,mex_name);

Param = checkInputs(Param,'REXstarJGG');

if Param.n_constraint == 0
%     Param.interimreportfun = @(x,f) x;
    Param.interimreportfun = @(x,f) interimreportfun_PE(x,f,model,mst,mex_name);
else
%     Param.interimreportfun = @(x,f) x;
    Param.interimreportfun = @(x,f,phi,g) interimreportfun_PE(x,f,phi,g,model,mst,mex_name);
end

% best = RCGA_PE_Main(Param,@JGG);
best = RCGA_Main(Param,@JGG);

x = Param.decodingfun(best.gene);

if isSBmodel(model)
    st_model = struct(model);
    for i = 1 : length(st_model.parameters)
        st_model.parameters(i).value = x(i);
    end
    optimizedmodel = SBmodel(st_model);
else
    optimizedmodel = [];
end

