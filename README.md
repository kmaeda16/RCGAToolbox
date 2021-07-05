# RCGAToolbox

Thank you for using RCGAToolbox!

RCGAToolbox is a MATLAB toolbox that contains two real-coded genetic algorithms (RCGAs): the unimodal normal distribution crossover with minimal generation gap (UNDX/MGG) and the real-coded ensemble crossover star with just generation gap (REXstar/JGG). The stochastic ranking method is implemented to efficiently handle constrained optimization problems. RCGAToolbox not only provides access to RCGAs but also several useful features for parameter estimation in systems biology. For demonstrations, run scripts in RCGAToolbox/doc/demo. For usage, type "help function_name" in the MATLAB Command Window. ***The complete user guide is provided as a Supplementary Material of the original paper (See Citation below)***.

## License

RCGAToolbox is distributed under GNU General Public License v3.0. For academic usage, RCGAToolbox is free. For other usages, please contact the author(s).

## Release Notes

- Apr 27 2021: RCGAToolbox-1.2: PEtab support. Friendlier error messages.
- Feb 13 2021: RCGAToolbox-1.1: Bug fix, Improved GUIs, Linux and macOS support.
- Dec  9 2020: RCGAToolbox-1.0.


## Requirements

- ***MATLAB R2016a or later***. We confirmed that RCGAToolbox runs on (i) Windows 10 (2004) with MATLAB R2016a, (ii) SUSE Linux Enterprise Server 11 (x86_64) with MATLAB R2016a, and (iii) macOS Big Sur (11.3.1, Intel CPU) with MATLAB R2021a.
- Optional requirements
    - ***Parallel Computing Toolbox*** is required for parallel computation (opts.n_par > 1). It is not required for sequential computation.
    - ***Optimization Toolbox*** is required for local optimization using fmincon (opts.local = 1). It is not required if the local optimization function is not used.
    - ***IQM Tools*** (formerly known as SBToolbox2/SBPD) are required for handling Systems Biology Markup Language (SBML) and a fast simulation (fast_flag = 2).
    - ***SundialsTB*** is required for a fast simulation with CVODE (fast_flag = 1).


## Installation

1. Download RCGAToolbox from https://github.com/kmaeda16/RCGAToolbox.
2. Run the installation script `RCGAToolbox_Install.m` under the directory `RCGAToolbox/install/`.

## Uninstallation

1. Run the uninstallation script `RCGAToolbox_Uninstall.m` under the directory `RCGAToolbox/install/`.
2. Delete the directory RCGAToolbox.

## Diagnosis

`RCGAToolbox/install/RCGAToolbox_Diagnosis.m` is the self-diagnosis script that checks whether the RCGAToolbox is properly installed. It also tests the RCGAToolbox functions that depend on optional toolboxes. For the diagnosis, run `RCGAToolbox_Diagnosis.m` under the directory `RCGAToolbox/install/`.

## Directories

- `doc`: User guide (not yet), demo scripts, and benchmark scripts are included. For details in the benchmark experiments, see the original article (Maeda et al., 2021).
- `install`: Installation, uninstallation, and diagnosis scripts are included.
- `source`: Source codes are included.

## Screenshot

<img src="Screenshot.png" height="320px"/>

## Citation

If you make use of RCGAToolbox in any publications, we kindly ask that the following paper is cited:

<a href="https://doi.org/10.1101/2021.02.15.431062">Maeda K, Boogerd FC, Kurata H: RCGAToolbox: A real-coded genetic algorithm software for parameter estimation of kinetic models. bioRxiv 2021:2021.2002.2015.431062.</a>


Any suggestions and bug reports are welcome. Please contact KM.


-------------------------------
Kazuhiro Maeda (Kyushu Institute of Technology)
