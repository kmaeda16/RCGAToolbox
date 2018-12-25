function [ x, optimizedmodel ] = REXstarJGG_PE(model,mst,fast_flag,Param,fitnessfun_PE)

if isSBmodel(model)
    if fast_flag == 1
        temp = SBstruct(model);
        mex_name = strcat(temp.name,'_mex');
        clear(mex_name);
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
        if fast_flag == 1
            temp = SBstruct(model);
            mex_name = strcat(temp.name,'_mex');
            clear(mex_name);
            SBPDmakeMEXmodel(model,mex_name);
        else
            mex_name = [];
        end
    end
end

if ~isSBmeasurement(mst)
    fprintf('Reading %s ...',mst);
    mst = SBmeasurement(mst);
    fprintf(' Finished.\n');
end


Param.fitnessfun = @(x) fitnessfun_PE(x,model,mst,mex_name);
interimreportfun = Param.interimreportfun;
Param.interimreportfun = @(elapsedTime,generation,Param,Population,best) interimreportfun(elapsedTime,generation,Param,Population,best,model,mst,mex_name,fast_flag);

Param = checkInputs(Param,'REXstarJGG');

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

