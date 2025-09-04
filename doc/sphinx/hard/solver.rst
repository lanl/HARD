.. |br| raw:: html

   <br />

.. _hard_solver_input:

HARD Solver Details
*******************

TBD

HARD Solver Parameters
**********************

All the options below are present under `linear_solver`:

`maxiter`
  Maximum number of iterations for the solver. Default is `50`.

`rtol`
  The residual tolerance for the solver. Default is `1e-12`.

`use_zero_guess`
  Boolean to define if the initial guess is 0. Default is `true`.

`flecsolve_coarse_grid`
  Boolean to define if we use the coarse grid in flecsolve. Default is `true`.

`jocabi_iterations`
  Define the maximum number of Jacobi iterations. Default is `100`.

`full_multigrid`
  Boolean to define if we use the full multigrid solver. Default is `true`.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
