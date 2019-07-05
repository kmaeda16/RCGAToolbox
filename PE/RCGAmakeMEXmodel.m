function mexfilename = RCGAmakeMEXmodel(model, mexfilename, doNOTcompileFlag)

if ~exist('doNOTcompileFlag','var')
    doNOTcompileFlag = 0;
end

if isIQMmodel(model)
    sbm = model;
    if ~exist('odefilename','var')
        st_sbm = struct(sbm);
        mexfilename = strcat(st_sbm.name,'_mex');
        mexfilename = regexprep(mexfilename,'\W','');
    end
elseif ischar(model)
    sbm = IQMmodel(model);
    if ~exist('odefilename','var')
        [ ~, filename, ~] = fileparts(model);
        mexfilename = strcat(filename,'_mex');
    end
else
    error('model must be an IQMmodel or file name!');
end


IQMmakeMEXmodel(sbm,mexfilename,1);

if doNOTcompileFlag == 1
    mexcompileIQM(mexfilename);
elseif doNOTcompileFlag ~= 0
    error('Unexpected doNOTcompileFlag!');
end
