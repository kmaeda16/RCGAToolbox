function RCGAcreateDecodingfun(PEtabParameterFile, DecodingfunFile)
% RCGAcreateDecodingfun creates a decoding function from  a PEtab parameter
% file.
% 
% [SYNTAX]
% RCGAcreateDecodingfun(PEtabParameterFile, DecodingfunFile)
% 
% [INPUT]
% PEtabParameterFile :  Name of PEtab parameter file. This is a TSV file
%                       with at least 6 fields: parameterId, 
%                       parameterScale, lowerBound, upperBound, 
%                       nominalValue, and estimate. For the format, see 
%                       Schmiester, PloS Comput Biol, 2021,
%                       DOI: 10.1371/journal.pcbi.1008646
% DecodingfunFile    :  Name of a decoding function file to be created from 
%                       the above PEtab parameter file.


%% Checking input arguments
if nargin ~= 2
    error('Incorrect number of input arguments.');
end


%% Preparation
[ ~, funcname, ~ ] = fileparts(DecodingfunFile);
filename = strcat(funcname,'.m');

T = tdfread(PEtabParameterFile);
temp = size(T.parameterId);
n_param = temp(1);


%% Opening Output File
out = fopen(filename,'w');
if out == -1
    warning('cannot open %s!\n',filename);
    return;
end


%% Function name
fprintf(out,'function param = %s(gene)\n', funcname);


%% Information
fprintf(out,'%% This function was created by RCGAToolbox RCGAcreateDecodingfun\n');
fprintf(out,'%% based on %s.\n',PEtabParameterFile);
fprintf(out,'%% Created: %s\n\n',date);


%% Mapping gene into param
for i = 1 : n_param;
    
    fprintf(out,'%% %s\n',strtrim(T.parameterId(i,:)));
    
    switch T.estimate(i)
        
        case 0
            fprintf(out,'param(%d) = %e;\n',i,T.nominalValue(i));
            
        case 1
            switch strtrim(T.parameterScale(i,:))
                case {'lin'}
                    fprintf(out,'lb = %e;\n',T.lowerBound(i));
                    fprintf(out,'ub = %e;\n',T.upperBound(i));
                    fprintf(out,'param(%d) = lb + ( ub - lb ) * gene(%d);\n',i,i);
                case {'log'}
                    fprintf(out,'lb = log(%e);\n',T.lowerBound(i));
                    fprintf(out,'ub = log(%e);\n',T.upperBound(i));
                    fprintf(out,'param(%d) = exp( lb + ( ub - lb ) * gene(%d) );\n',i,i);
                case {'log10'}
                    fprintf(out,'lb = log10(%e);\n',T.lowerBound(i));
                    fprintf(out,'ub = log10(%e);\n',T.upperBound(i));
                    fprintf(out,'param(%d) = 10 ^ ( lb + ( ub - lb ) * gene(%d) );\n',i,i);
                otherwise
                    error('Unexpected parameterScale in %s. lin, log, and log10 are supported.',PEtabParameterFile);
            end
            
        otherwise
            error('Unexpected estiamte in %s. It should be 0 (fixed) or 1 (subject to estimation).',PEtabParameterFile);
    end
    
    fprintf(out,'\n');
    
end


%% Print Path and File name
fprintf('Decodingfun "%s" was created in "%s".\n',filename,pwd);


%% Closing Output File
fclose(out);
