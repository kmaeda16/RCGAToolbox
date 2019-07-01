function f = SSR(Simulation, x, modelfun, mst, opts)


%% Setting experimental data
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


%% Running simulation
y0 = modelfun();
[ ~, X ] = feval(Simulation, modelfun, tspan, y0, x', opts);


%% If simulation faild, return f = 1e+10
if max(max(isnan(X)))
    warning('Simulation failed!');
    f = 1e+10;
    return;
end


%% Preparing x_sim and x_exp
x_sim = X;
[n_row, n_col] = size(x_sim);
x_exp = zeros([n_row n_col]);
for i = 1 : n_col
    x_exp(:,i) = mst.data(i).values;
end


%% Calculating SSR
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

