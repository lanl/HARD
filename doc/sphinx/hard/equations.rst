.. |br| raw:: html

   <br />

.. _hard_equations:


Basic Equations
***************

We follow the comoving frame formulation of the radiation hydrodynamics equations in the diffusion limit, which can be written in Einstein notation in the form

.. math::
   \partial_t \mathbf{U} + \partial_j \mathbf{F}^j = \mathbf{S}

with the evolved variables

.. math::
   \mathbf{U} = \left[ \rho, \rho v^i, e + \frac{1}{2} \rho v^2, E \right]^T,

fluxes

.. math::
   \mathbf{F}^j = \left[ \begin{matrix}
        \rho v^j \\
        \rho v^i v^j + p \delta_{ij} \\
        \left(e + \frac{1}{2} \rho v^2 + p\right) v^j \\
        E v^j \\
        \end{matrix} \right],

and source terms

.. math::
   \mathbf{S} = \left[ \begin{matrix}
   0 \\
   f^i \\
   f^i v_i + \dot{q} \\
   - \partial_i F_r^i - P_r^{ij} \nabla_i v_j - \dot{q} \\
   \end{matrix} \right],

where :math:`\rho` is the fluid mass density, :math:`v^i` is the component *i* of the fluid velocity vector, :math:`e` is the fluid internal energy density, :math:`p` is the fluid pressure, :math:`E` is radiation energy density, :math:`F_r^i` is the *i*-th component of the radiation energy flux vector, :math:`P_r^{ij}` is the *i,j* element of the radiation pressure tensor, :math:`f^i` is the *i*-th component of the radiation force density vector

.. math::
   f^i = \frac{\kappa \rho}{c} F_r^i ,

and :math:`\dot{q}` is the radiative heating rate

.. math::
   \dot{q} = \frac{dq}{dt} = c \kappa \rho (E - a T^4) ,

where :math:`\kappa` is mean opacity, :math:`c` is the speed of light, :math:`a` is the radiation constant, and :math:`T` is the fluid temperature.

Flux Limited Diffusion
~~~~~~~~~~~~~~~~~~~~~~~~

For closure of the evolution system for the fluid, we need to specify the functional form of the fluid pressure :math:`p`. Here we adopt the singularity-EOS, allowing us to use any material in evolution.

The radiative transfer equation — the last component of the evolution equation :math:`\partial_t \mathbf{U} + \partial_j \mathbf{F}^j = \mathbf{S}` — requires additional closures for radiation energy flux density :math:`F^i` and radiation pressure :math:`P_r^{ij}`.

In the optically thick diffusion limit, they are given by:

.. math::
   F^i = - \frac{c}{3\kappa\rho} \nabla E,

.. math::
   P_r^{ij} = \frac{1}{3} E \delta_{ij},

provided the mean free path of photons is much shorter than the length scale of interest.

Evolution equations are valid in the diffusion (optically thick) limit. However, it is possible to introduce a bridge law so that the radiation energy flux density :math:`F_r^i` recovers the correct norm in optically thin regions.

.. math::
   F_r^i = - \lambda \frac{c}{\kappa \rho} \nabla E

:math:`\lambda` is a radiation energy flux limiter function:

.. math::
   \lambda = \frac{2 + R}{6 + 3R + R^2}

.. math::
   R \equiv \frac{|\nabla E|}{\rho \kappa E}

.. math::
   \frac{P^{ij}}{E} = \left(\frac{1 - f}{2}\right) \delta_{ij}
      + \left(\frac{3f - 1}{2}\right) n^i n^j

where :math:`n^i \equiv \nabla E / |\nabla E|` and the quantity :math:`f` is defined as:

.. math::
   f = \lambda + \lambda^2 R^2

In the optically thick limit, :math:`R \to 0`, :math:`\lambda \to 1/3`, :math:`f \to 1/3`, and the Eddington approximation :math:`P^{ij} = (E/3)\delta_{ij}` is recovered.

In the optically thin limit, the photon field streams toward the direction of the local radiation energy density gradient at the speed of light.

These prescriptions are not a physically exact implementation of realistic light-matter interactions but serve as a computationally efficient proxy for bridging optically thin and thick limits.

As a result, the radiative transfer equation becomes:

.. math::
   \partial_t E + \partial_j (Ev^j) =
        \partial_i(D\, \partial_i E) - P_r^{ij} \nabla_i v_j - \dot{q},

with :math:`D = c\lambda / (\kappa\rho)`.

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :


