# RCGAToolbox

Thank you for using RCGAToolbox!

RCGAToolbox is a MATLAB toolbox for parameter estimation in systems biology, which contains two real-coded genetic algorithms (RCGAs), viz. the Unimodal Normal Distribution Crossover with Minimal Generation Gap (UNDX/MGG) and the Real-coded Ensemble Crossover star with Just Generation Gap (REXstar/JGG). The stochastic ranking method is implemented to handle constrained optimization problems efficiently. RCGAToobox provides not only the access to RCGAs but also to several useful features for parameter estimation. For demonstrations, run scripts in RCGAToolbox/doc/demo. For usage, type "help function_name" in the MATLAB Command Window. ***The complete user guide is provided as a Supplementary Material of the original paper (See Citation below)***.

## License

RCGAToolbox is distributed under GNU General Public License v3.0. For academic usage, RCGAToolbox is free. For other usages, please contact the author(s).

## Release Note

Dec 9 2020: RCGAToolbox-1.0 released.

## Requirements

- MATLAB R2016a or later [We confirmed that RCGAToolbox runs on Windows 10 (1909) and SUSE Linux Enterprise Server 11 (x86_64)].
- Optional requirements
    - Parallel Computing Toolbox is required for parallel computation (opts.n_par > 1). For serial computation, Parallel Computing Toolbox is not required.
    - Optimization Toolbox is required for local optimization using fmincon (opts.local = 1). It is not required if you do not use the local optimization function.
    - IQM Tools (formerly known as SBToolbox2/SBPD/SBPOP) is required for handling SBML (Systems Biology Markup Language) and a fast simulation (fast_flag = 2).
    - SundialsTB is required for a fast simulation with CVODE (fast_flag = 1).


## Installation

1. Download RCGAToolbox from https://github.com/kmaeda16/RCGAToolbox.
2. Run the installation script RCGAToolbox_Install.m under the directory RCGAToolbox/install/.

## Uninstallation

1. Run the uninstallation script RCGAToolbox_Uninstall.m under the directory RCGAToolbox/install/.
2. Delete the directory RCGAToolbox.

## Diagnosis

RCGAToolbox/install/RCGAToolbox_Diagnosis.m is the diagnosis script which checks whether RCGAToolbox is properly installed. It also tests the RCGAToolbox functions that depend on optional toolboxes. For the diagnosis, run RCGAToolbox_Diagnosis.m under the directory RCGAToolbox/install/.

## Citation

Maeda and Kurata, RCGAToolbox: a toolbox for real-coded genetic algorithms for parameter estimation of kinetic models


Good luck! Any suggestions and bug reports are welcome. Please contact KM.


-------------------------------
Kazuhiro Maeda,
Kyushu Institute of Technology
