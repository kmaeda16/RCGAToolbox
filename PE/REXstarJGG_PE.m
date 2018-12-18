function [ optimizedmodel, x ] = REXstarJGG_PE(model,measurment,opts)


if isSBmodel(model)
    if opts.fast == 1
        temp = SBstruct(model);
        SBPDmakeMEXmodel(model,strcat(temp.Name,'_mex'));
    end
else
    switch exist(model,'file')
        case 3
            opts.fast = 1;
        case 2
            model = SBmodel(model);
            if opts.fast == 1
                temp = SBstruct(model);
                SBPDmakeMEXmodel(model,strcat(temp.Name,'_mex'));
            end
        otherwise
            error('model should be an SBmodel, SBML file, or MEX file');
    end
end

opts.model = model;
opts.mst = measurment;
opts.ub = opts.ub;
opts.lb = opts.lb;
opts.mex_name = 'testtest';

best = RCGA_PE_Main(opts,@JGG);

x = best.gene .* ( opts.ub - opts.lb ) + opts.lb;

st_model = struct(model);
for i = 1 : length(st_model.parameters)
    st_model.parameters(i).value = x(i);
end

optimizedmodel = SBmodel(st_model);
