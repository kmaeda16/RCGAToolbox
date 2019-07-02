function f = SSR(param, Simulation, modelfun, mst, simopts)
% SSR calculates the sum of squared resudials between simulation and
% experimental data
% 
% [SYNTAX]
% f = SSR(x, Simulation, modelfun, mst)
% f = SSR(x, Simulation, modelfun, mst, simopts)
% 
% [INPUT]
% param     :  Parameter value vector
% Simulation:  Function handle for Simulation_*
% modelfun  :  Function handle for model (odefun or mex)
% mst       :  Experimental data (IQMmeasurement)
% simopts   :  Structure with integrator options. Fields depend on
%              Simulation_*. See 'help Simulation_'.
% 
% [OUTPUT]
% f:  Objective function value (scaler)


%% Handling inputs
if nargin == 4
    simopts = [];
end


%% Converting IQMmeasuremen into structure
mst = struct(mst);


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
[ ~, Y ] = feval(Simulation, modelfun, tspan, y0, param, simopts);


%% If simulation faild, return f = 1e+10
if max(max(isnan(Y)))
    warning('Simulation failed!');
    f = 1e+10;
    return;
end


%% Preparing Y_sim and Y_exp
Y_sim = Y;
[n_row, n_col] = size(Y_sim);
Y_exp = zeros([n_row n_col]);
for i = 1 : n_col
    Y_exp(:,i) = mst.data(i).values;
end


%% Calculating SSR
f = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(Y_exp(i,j))
            % f = f + abs( Y_sim(i,j) - Y_exp(i,j) );
            % f = f + abs( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) );
            f = f + ( Y_sim(i,j) - Y_exp(i,j) ) ^ 2;
            % f = f + ( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) ) ^ 2;
        end
    end
end
