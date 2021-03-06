function [ fout_name, flag ] = RCGAreplaceWords(fin_name, fout_name, list)
% RCGAreplaceWords replaces "word" in list with "word0", which avoids
% errors caused by wrongly used reserved words. Use this function when the
% function mexcompileIQM gives you errors because of C-language-reserved
% words wrongly used in a C-source code.
% 
% [SYNTAX]
% [ fout_name, flag ] = RCGAreplaceWords(fin_name)
% [ fout_name, flag ] = RCGAreplaceWords(fin_name, fout_name)
% [ fout_name, flag ] = RCGAreplaceWords(fin_name, fout_name, list)
% 
% [INPUT]
% fin_name  :  Name of the input file.
% fout_name :  Name of the output file (default: 'inputfile_mod.txt' if 
%              fin_name = 'inputfile.txt'.
% list      :  Cell array of words to be replaced (default: Some of 
%              reserved words for C language. See the codes below).
% 
% [OUTPUT]
% fout_name :  Name of the output file.
% flag      :  flag = 0 for successs, and flag = 1 for error. 


%% Set default input variables
if ~exist('fout_name','var')
    [pathstr, name, ext] = fileparts(fin_name);
    fout_name = strcat(pathstr,name,'_mod',ext);
end

if ~exist('list','var')
    list = {'default','auto','const','do','enum','extern','goto', ...
    'register','signed','sizeof','short','typedef','unsigned', ...
    'volatile'};
end


%% Opening files
fin = fopen(fin_name,'r');
if fin == -1
    warning('cannot open %s!\n',fin_name);
    flag  = 1;
    return;
end

fout = fopen(fout_name,'w');
if fout == -1
    warning('cannot open %s!\n',fout);
    flag  = 1;
    return;
end


%% Replacing words in list
while ~feof(fin)
    line = fgetl(fin);
    for i = 1 : length(list)
        expr = sprintf('%s%s%s','(^|\W)',char(list(i)),'($|\W)');
        replace = sprintf('$1%s%s$2',char(list(i)),'0');
        line = regexprep(line,expr,replace);
    end
    fprintf(fout,'%s\n',line);
end


%% Finishing 
fclose(fin);
fclose(fout);
flag = 0;

