.. |br| raw:: html

   <br />

.. _hard_problem:

Example Problem
***************

The application codes in this example solve the *Euler* equations of
compressible gas dynamics. In one spatial dimension, the Euler equations
are given by

.. math::

   \left[
     \begin{gathered}
       \rho \\
       \rho u \\
       E
     \end{gathered}
   \right]_t
   +
   \left[
     \begin{gathered}
       \rho u \\
       \rho u^2 + p \\
       \left(E + p\right)u
     \end{gathered}
   \right]_x
   = 0,

where :math:`\rho` is the density, :math:`u` is the velocity, :math:`E` is the
total engery, and :math:`p` is the pressure. This is generally referred to as
the *conservation* form, as it emphasises the property that the Euler equations
are mass, momentum, and total energy conserving.

In three dimensions, the conservation form is

.. math::

   q_t + f(q)_x + g(q)_y + h(q)_w = 0,

where

.. math::

   q =
   \left[
     \begin{gathered}
       \rho \\
       \rho u \\
       \rho v \\
       \rho w \\
       E
     \end{gathered}
   \right],
   f(q) =
   \left[
     \begin{gathered}
       \rho u \\
       \rho u^2 + p \\
       \rho uv \\
       \rho uw \\
       \left(E+p\right)u
   \end{gathered}
   \right],
   g(q) =
   \left[
     \begin{gathered}
       \rho v \\
       \rho uv \\
       \rho v^2 + p \\
       \rho vw \\
       \left(E+p\right)v
     \end{gathered}
   \right],
   h(q) =
   \left[
     \begin{gathered}
       \rho w \\
       \rho uw \\
       \rho vw \\
       \rho w^2 + p \\
       \left(E+p\right)w
     \end{gathered}
   \right].

MUSCL Scheme
~~~~~~~~~~~~

The MUSCL scheme (*Monotonic Upstream-centered Scheme for Conservation
Laws*) is a second-order (in space) finite volume method that is *total
variation diminishing (TVD)* and can provide accurate simulations of
shock dynamics. A succinct, semi-discrete representation of the scheme
can be written

.. math::

   \frac{du_i}{dt} + \frac{1}{\Delta x_i}
   \left[
     F^{*}_{i+1/2} - F^{*}_{i-1/2}
   \right],\label{semi}\tag{1}

where :math:`F^{*}_{i\pm 1/2}` are numerical fluxes that are a
non-linear mix of first and second-order slope reconstructions depending
on the gradients close to the current cell, and can be written in the
form

.. math::

   F^{*}_{i\pm 1/2} = f^{low}_{i\pm 1/2} -
   \phi(r_i)
   \left(
     f^{low}_{i\pm 1/2} - f^{high}_{i\pm 1/2}
   \right)

where :math:`r_i` is the ratio of successive gradient approximations and
:math:`\phi: r_i \rightarrow [0,1]` is a limiter function that selects a
realistic spatial derivative value. One example of such a limiter
function is the *minmod* limiter

.. math::

   \phi_{mm}(r_i) = max[0, min(1,r_i)];
   \;\;\lim_{r_i \rightarrow \infty}\phi_{mm}(r_i) = 1.

In smooth areas of the flow, :math:`\phi_{mm}(r_i)` will be close to
one, which will select the higher-order approximation. Near a shock,
:math:`r_i` will approach infinity and :math:`\phi_{mm}(r_i)` will be
close to zero and the low-order approximation will be selected.

Hancock Time Evolution
~~~~~~~~~~~~~~~~~~~~~~

To fully discretize (:math:`\ref{semi}`), we can use the Hancock method,
which is second-order accurate in time. Hancock is a predictor-corrector
method that uses the averages of the reconstructed input quantities to
approximate the flux for :math:`t+\Delta t/2` in the predictor
step. These intermediate quantities are then reconstructed to the faces
again, but a *Riemann* solver is used to compute the corrector step.

.. note::

   The MUSCL-Hancock method documented here is second-order accurate in
   space and time.

Quantity reconstruction (extrapolation) to faces uses the rule

.. math::

   q^n_{i\pm 1/2} = q^n_i \pm
     \frac{\Delta x_i}{2}
     \left(\frac{\partial q^n_i}{\partial s}\right),\label{recon}\tag{2}

where :math:`\partial q^n_i/\partial s` is the flux limited slope.

The predictor step computes intermediate quantities :math:`q^{*}` using

.. math::

   q^{*}_i = q^n_i - \frac{\Delta t}{2\Delta x_i}
   \left(q^n_{i-1/2} - q^n_{i+1/2}\right).

These quantities are then reconstructed to the faces again by applying
(:math:`\ref{recon}`). The corrector step is

.. math::

   q^{n+1}_i = q^n_i + \frac{\Delta t}{\Delta x}
   \left( F_{i-1/2} - F_{i+1/2} \right),

where :math:`F()` is a Riemann solver. The code in this specialization
example uses the HLL approximate Riemann sovler.

.. math::

   F^{hll} =
   \frac
   {S_R F_L - S_L F_R + S_R S_L\left(U_R - U_L\right)}
   {S_R - S_L}

where :math:`S_T` and :math:`S_H` are the min and max wave speeds taken
at the tail and head of the cell interfaces respectively, :math:`U_L`
and :math:`U_R` are the respective intermediate quantities, and
:math:`F_L` and :math:`F_R` are the flux functions evaluated at the
respective intermediate quantities.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
