To start benchmark experiments, 
1. For Windows, modify mex_threestep.m and mex_hiv.m, and rum them to compile threestep_con_mex and hiv_con_mex. For Linux and macOS, threestep_con_mex and hiv_con_mex are provided. So, you don't need to compile them.
2. Run run_Benchmark.m in the directories MATLAB_GA, eSS, and RCGA.

- For the function batch, Parallel Computing Toolbox is required.
- For the function ga, Global Optimization Toolbox is required.
- For the function ess_kernel, MEIGO is required.
- For compilation of threestep_con_mex and hiv_con_mex, sundials-2.4.0 is required.
