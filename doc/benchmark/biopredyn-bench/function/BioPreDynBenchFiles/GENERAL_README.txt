This is the GENERAL README file for the BioPreDyn-Bench benchmarks.
http://www.iim.csic.es/~gingproc/biopredynbench/

Created:	17/06/2014
Last modified:	20/06/2014
Contact:	gingproc@iim.csic.es


IMPLEMENTATIONS

We provide the following formats:

(1) ready-to-run implementations for Matlab, AMIGO (Matlab toolbox) and COPASI
(2) core implementations in C (which only need the user to write a small program as 
optimization driver, therefore allowing the use of any optimization code which can be 
interfaced with C)
(3) partial implementations (systems dynamics plus initial conditions) in SBML 
(so they can be imported into any of the many software packages that support this format)

In the last case, to fully define the parameter estimation problems, users would 
also need to align the data with the observed states and build the cost functions to be 
optimized.

TABLE OF COMPATIBILITY

We indicate below, for each benchmark implementation, the OPERATING SYSTEM 
where it can be executed and ADDITIONAL SOFTWARE REQUIREMENTS (indicated 
inside parentheses).

The following shortnames are used in the table:

Win: Windows; tested with Windows 7 64-Bit
Lin: Linux; tested with openSUSE 64-Bit
OSX: Mac OSX
ML32: Matlab 32-Bit version, tested with Matlab R2011 32-Bit for Windows
ML64: Matlab 64-Bit version, tested with Matlab R2008a 64-Bit for Linux
MCR:  Matlab Compiler Runtime, available at http://www.mathworks.es/products/compiler/mcr/


+----------------+----------------------+----------------------+----------------------+
| Implementation |          B1          |          B2          |          B3          |
+----------------+----------------------+----------------------+----------------------+
| AMIGO          | Win(ML32), Lin(ML64) | Win(ML32), Lin(ML64) | Win(ML32), Lin(ML64) |
| C              | Win(MCR),Lin(MCR)    | Win(MCR),Lin(MCR)    | Win(MCR),Lin(MCR)    |
| COPASI         | Win,Lin,OSX          | Win,Lin,OSX          | Win,Lin,OSX          |
| Matlab         | Win(ML32), Lin(ML64) | Win(ML32), Lin(ML64) | Win(ML32), Lin(ML64) |
| SBML           | Win,Lin,OSX          | Win,Lin,OSX          | Win,Lin,OSX          |
+----------------+----------------------+----------------------+----------------------+


+----------------+----------------------+----------------------+-----------------+
| Implementation |          B4          |          B5          |       B6        |
+----------------+----------------------+----------------------+-----------------+
| AMIGO          | Win(ML32), Lin(ML64) | Win(ML32), Lin(ML64) | Lin(ML32,ML64)  |
| C              | Win(MCR),Lin(MCR)    | Win(MCR),Lin(MCR)    | Lin(MCR)        |
| COPASI         | Win,Lin,OSX          |                      |                 |
| Matlab         | Win(ML32),Lin(ML64)  | Win(ML32),Lin(ML64)  | Lin(ML32,ML64)  |
| SBML           | Win,Lin,OSX          | Win,Lin,OSX          |                 |
+----------------+----------------------+----------------------+-----------------+

Additional notes:
- Matlab Optimization Toolbox might be needed for certain local solvers
- C implementations: see further details in the corresponding readme files
  regarding further requirements, installation and compilation instructions
