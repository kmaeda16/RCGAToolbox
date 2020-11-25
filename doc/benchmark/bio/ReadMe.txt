To start benchmark experiments, 
1. Compile threestep_con_mex and hiv_con_mex by running mex_threestep.m and mex_hiv.m with modification.
2. Run run_Benchmark.m in the directories MATLAB_GA, eSS, and RCGA.

- For the function batch, Parallel Computing Toolbox is required.
- For the function ga, Global Optimization Toolbox is required.
- For the function ess_kernel, SSm GO Toolbox is required.
- For compilation of threestep_con_mex and hiv_con_mex, sundials-2.4.0 is required.
