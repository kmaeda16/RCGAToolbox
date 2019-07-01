function [f, g]= SSR(Simulation, x, model, mst, opts)
% function f = SSR_sbml(x, mex_name, mst, opts)

if length(mst) > 2
    warning('Measurement has multiple measurment deta sets, but only the first data set will be used for fitness calculation.');
end
mst = struct(mst{1});

%% Checking time errors
t0 = 0;
if t0 < mst.time(1)
    tspan = [ t0  mst.time ];
elseif t0 == mst.time(1)
    tspan = mst.time;
else
    error('Time of the first experimental datapoint should be AFTER or EQUAL TO time 0');
end

%%
y0 = model();
[ ~, X ] = feval(Simulation, model, tspan, y0, x', opts);

if max(max(isnan(X)))
    f = 1e+10;
    return;
end

x_sim = X;
x_exp = zeros(size(x_sim));
for i = 1 : length(mst.data)
    x_exp(:,i) = mst.data(i).values;
end


[n_row, n_col] = size(x_exp);
f = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(x_exp(i,j))
%             f = f + abs( x_sim(i,j) - x_exp(i,j) );
%             f = f + abs( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) );
            f = f + ( x_sim(i,j) - x_exp(i,j) ) ^ 2;
%             f = f + ( ( x_sim(i,j) - x_exp(i,j) ) / x_exp(i,j) ) ^ 2;
        end
    end
end
g = 0;