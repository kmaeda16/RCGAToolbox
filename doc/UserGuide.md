# Function Lists

- [odestb](##odestb)

---

## odestb

### Description

odestb is a wrapper function that enables to use CVode provided by
SundialsTB in the same way as MATLAB ODE solvers.

### Syntax

- [ T, Y ] = odestb(odefun, tspan, y0)
- [ T, Y ] = odestb(odefun, tspan, y0, options)

### Input

- odefun :  ODEFUN file.
- tspan  :  [t0, tf] or [t0, t1, ..., tf].
- y0     :  Initial value vector.
- options:  Structure with integrator options.
  - options.LMM: Linear Multistep Method (default: 'BDF')
  - options.NonlinearSolver: Type of nonlinear solver used (default: 'Newton')
  - options.AbsTol:  Absolute tolerance (default: 1e-6)
  - options.RelTol:  Relative tolerance (default: 1e-4)
  - options.MinStep: Minimum stepsize (default: 0)
  - options.MaxStep: Maximum stepsize (default: inf)
  - options.MaxNumSteps: Maximum number of steps (default: 500)
    For other fields, see SundialsTB documentation.

### Output

- T      :  Column vector of timepoints.
- Y      :  Variable matrix. Each column corresponds to each variable. Each row of Y corresponds to each row of T.

---

## odestx

### Description

odestb is a wrapper function that enables to use CVode provided by
SundialsTB in the same way as MATLAB ODE solvers.

### Syntax

- [ T, Y ] = odestb(odefun, tspan, y0)
- [ T, Y ] = odestb(odefun, tspan, y0, options)

### Input

- odefun :  ODEFUN file.
- tspan  :  [t0, tf] or [t0, t1, ..., tf].
- y0     :  Initial value vector.
- options:  Structure with integrator options.
  - options.LMM: Linear Multistep Method (default: 'BDF')
  - options.NonlinearSolver: Type of nonlinear solver used (default: 'Newton')
  - options.AbsTol:  Absolute tolerance (default: 1e-6)
  - options.RelTol:  Relative tolerance (default: 1e-4)
  - options.MinStep: Minimum stepsize (default: 0)
  - options.MaxStep: Maximum stepsize (default: inf)
  - options.MaxNumSteps: Maximum number of steps (default: 500)
    For other fields, see SundialsTB documentation.

### Output

- T      :  Column vector of timepoints.
- Y      :  Variable matrix. Each column corresponds to each variable. Each row of Y corresponds to each row of T.