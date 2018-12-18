function f = mySSR(model,mex_name,mes)
% function f = mySSR(model,mes,w)

% model = SBmodel('SBMLexampleLevel2.xml');
% mes = SBmeasurement('MeasurementExample.xls');

mes = struct(mes{1});

% if isSBmodel(model)
%     output = SBsimulate(model,mes.time);
% else
    temp = struct(model);
    for i = 1 : length(temp.parameters)
        parameters{i} = temp.parameters(i).name;
        parametervalues(i) = temp.parameters(i).value;
    end
    parametervector = makeparamvecSBPD(mex_name,parameters,parametervalues);
%     output = SBPDsimulate(mex_name,mes.time,[],[],parametervector);
%     output = SBPDsimulate(mex_name,mes.time,[],[],zeros(1,9));
%     output = testtest(mes.time,[],parametervector);
try
    output = feval(mex_name,mes.time,[],parametervector);
catch
    f = 1e+10;
    return;
end
% end

x_sim = output.statevalues;
w = ones(size(x_sim));

if max(max(isnan(x_sim)))
    f = 1e+10;
    return;
end

x_exp = zeros(size(x_sim));
for i = 1 : length(mes.data)
    x_exp(:,i) = mes.data(i).values;
end

% f = sum( sum( ( ( x_sim - x_exp ) ./ x_exp ) .^ 2 ) );
% f = sum( sum( ( x_sim - x_exp ) .^ 2 ) );
% g = 0;

[n_row, n_col] = size(x_exp);
f = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(x_exp(i,j))
            f = f + w(i,j) * abs( x_sim(i,j) - x_exp(i,j) );
%             f = f + w(i,j) * abs( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) );
%             f = f + w(i,j) * ( x_sim(i,j) - x_exp(i,j) ) ^ 2;
%             f = f + w(i,j) * ( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) ) ^ 2;
        end
    end
end

% plot(mes.time,x_exp,'o');
% hold on;
% ax = gca;
% ax.ColorOrderIndex = 1;
% plot(output.time,x_sim);
% hold off;

