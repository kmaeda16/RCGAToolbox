function [f, g] = wrapper_threestep_con_mex(x)

x = 10 .^ x;

[f, g] = threestep_con_mex(x);
