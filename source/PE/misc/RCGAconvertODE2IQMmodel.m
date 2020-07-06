function model = RCGAconvertODE2IQMmodel(ODEfileName)


% fileID = fopen('ExampleModel.m','r');
fileID = fopen(ODEfileName,'r');

ModelStruct = struct(IQMmodel);

tline = fgetl(fileID);
while ischar(tline)
    
    % === BEGIN NAME ===
    if strfind(tline,'=== BEGIN NAME ===')
        tline = fgetl(fileID);
        while ischar(tline)
            if strfind(tline,'=== END NAME ===')
                break;
            end
            tline = regexprep(tline,'^%\s*','');
            ModelStruct.name = strtrim(tline);
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN NOTES ===
    if strfind(tline,'=== BEGIN NOTES ===')
        tline = fgetl(fileID);
        while ischar(tline)
            if strfind(tline,'=== END NOTES ===')
                break;
            end
            tline = regexprep(tline,'^%','');
            ModelStruct.notes = strcat(ModelStruct.notes,tline);
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN INITIAL CONDITION ===
    if strfind(tline,'=== BEGIN INITIAL CONDITION ===')
        tline = fgetl(fileID);
        index = 0;
        while ischar(tline)
            if strfind(tline,'=== END INITIAL CONDITION ===')
                break;
            end
            name = regexprep(tline,'\s*(.+)_0\s*=.+','$1');
            initialCondition = regexprep(tline,'.+=\s*(.+)\s*;','$1');
            index = index + 1;
            ModelStruct.states(index).name = strtrim(name);
            ModelStruct.states(index).initialCondition = str2double(initialCondition);
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN PARAMETERS ===
    if strfind(tline,'=== BEGIN PARAMETERS ===')
        tline = fgetl(fileID);
        index = 0;
        while ischar(tline)
            if strfind(tline,'=== END PARAMETERS ===')
                break;
            end
            name = regexprep(tline,'\s*(.+)\s*=.+','$1');
            value = regexprep(tline,'.+=\s*(.+)\s*;','$1');
            index = index + 1;
            ModelStruct.parameters(index).name = strtrim(name);
            ModelStruct.parameters(index).value = str2double(value);
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN VARIABLES ===
    if strfind(tline,'=== BEGIN VARIABLES ===')
        tline = fgetl(fileID);
        index = 0;
        while ischar(tline)
            if strfind(tline,'=== END VARIABLES ===')
                break;
            end
            name = regexprep(tline,'\s*(.+)\s*=.+','$1');
            formula = regexprep(tline,'.+=\s*(.+)\s*;','$1');
            index = index + 1;
            ModelStruct.variables(index).name = strtrim(name);
            ModelStruct.variables(index).formula = formula;
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN REACTIONS ===
    if strfind(tline,'=== BEGIN REACTIONS ===')
        tline = fgetl(fileID);
        index = 0;
        while ischar(tline)
            if strfind(tline,'=== END REACTIONS ===')
                break;
            end
            name = regexprep(tline,'\s*(.+)\s*=.+','$1');
            formula = regexprep(tline,'.+=\s*(.+)\s*;','$1');
            index = index + 1;
            ModelStruct.reactions(index).name = strtrim(name);
            ModelStruct.reactions(index).formula = formula;
            tline = fgetl(fileID);
        end
    end
    
    % === BEGIN BALANCE ===
    if strfind(tline,'=== BEGIN BALANCE ===')
        tline = fgetl(fileID);
        index = 0;
        while ischar(tline)
            if strfind(tline,'=== END BALANCE ===')
                break;
            end
            name = regexprep(tline,'\s*(.+)_dot\s*=.+','$1');
            formula = regexprep(tline,'.+=\s*(.+)\s*;','$1');
            index = index + 1;
            % ModelStruct.states(index).name = strtrim(name);
            ModelStruct.states(index).ODE = formula;
            tline = fgetl(fileID);
        end
    end
    
    
    tline = fgetl(fileID);
end

fclose(fileID);


model = IQMmodel(ModelStruct);
