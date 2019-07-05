function odefilename = RCGAmakeODEmodel(model, odefilename)

if isIQMmodel(model)
    sbm = model;
    if ~exist('odefilename','var')
        st_sbm = struct(sbm);
        odefilename = strcat(st_sbm.name,'_odefun');
        odefilename = regexprep(odefilename,'\W','');
    end
elseif ischar(model)
    sbm = IQMmodel(model);
    if ~exist('odefilename','var')
        [ ~, filename, ~] = fileparts(model);
        odefilename = strcat(filename,'_odefun');
    end
else
    error('model must be an IQMmodel or file name!');
end

IQMcreateODEfile(sbm,odefilename);
