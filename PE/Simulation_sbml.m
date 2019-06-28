function  [ T, X ] = Simulation_sbml(x, mex_name, tspan, fast_flag, opts)

switch fast_flag
    case 0
        [ T, X ] = Simulation_odexx(x, mex_name, tspan, opts);
    case 1
        [ T, X ] = Simulation_stb(x, mex_name, tspan, opts);
    case 2
        [ T, X ] = Simulation_mex(x, mex_name, tspan, opts);
    otherwise
        error('Unexpected fast_flag!');
end
