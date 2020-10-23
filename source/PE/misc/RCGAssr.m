function ssr = RCGAssr(param, Simulation, model, mst, simopts)
% RCGAssr calculates the sum of squared resudials (SSR) between simulation
% and experimental data.
% 
% [SYNTAX]
% ssr = RCGAssr(param, Simulation, model, mst)
% ssr = RCGAssr(param, Simulation, model, mst, simopts)
% 
% [INPUT]
% param      :  Parameter value vector.
% Simulation :  Function handle for RCGAsimulateODEXX, RCGAsimulateSTB, or
%               RCGAsimulateMEX.
% model      :  Function handle for an ODE function (IQM Tools format) or a 
%               MEXed model.
% mst        :  Experimental data (An IQMmeasurement object).
% simopts    :  Solver option structure. The fields depend on fast_flag. 
%               For fast_flag = 0, 1, and 2, see 'help RCGAsimulateODEXX', 
%               'help RCGAsimulateSTB', 'help RCGAsimulateMEX', 
%               respectively.
% 
% [OUTPUT]
% ssr        :  Sum of squared resudials (SSR).


%% Handling inputs
if nargin == 4
    simopts = [];
end


%% Converting IQMmeasurement into structure
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
y0 = model();
[ ~, Y ] = feval(Simulation, model, tspan, y0, param, simopts);


%% If simulation faild, return ssr = 1e+10
[n_row, ~] = size(Y);
if max(max(isnan(Y))) || length(tspan) ~= n_row
    warning('Simulation failed!');
    ssr = 1e+10;
    return;
end


%% Preparing Y_sim and Y_exp
Y_sim = real(Y);
[n_row, n_col] = size(Y_sim);
Y_exp = zeros([n_row n_col]);
for i = 1 : n_col
    Y_exp(:,i) = mst.data(i).values;
end


%% Calculating SSR
ssr = 0;
for i = 1 : n_row
    for j = 1 : n_col
        if ~isnan(Y_exp(i,j))
            % ssr = ssr + abs( Y_sim(i,j) - Y_exp(i,j) );
            % ssr = ssr + abs( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) );
            ssr = ssr + ( Y_sim(i,j) - Y_exp(i,j) ) ^ 2;
            % ssr = ssr + ( ( Y_sim(i,j) - Y_exp(i,j) ) / Y_exp(i,j) ) ^ 2;
        end
    end
end
