.. |br| raw:: html

   <br />

.. _hard_test_problems:

Test Problems
************************************

We describe several test problems that users can test with HARD

Sod Shock Tube
~~~~~~~~~~~~~~~~~

The Sod shock tube is a standard test with a classical Riemann problem with the following initial parameters:

.. math::

   (\rho, v, p)_{t=0} =
   \begin{cases}
   (1.0, 0.0, 1.0) & \text{if} \;\; 0.0 < x \leq 0.5 \\
   (0.125, 0.0, 0.1) & \text{if} \;\; 0.5 < x < 1.0.
   \end{cases}

This leads to the development of a shock front, which propagates from high-density into low-density regions, and is followed by a contact discontinuity. A density rarefaction wave propagates into the high-density region.

Heating-Cooling
~~~~~~~~~~~~~~~~~



Temperature Induced Shock
~~~~~~~~~~~~~~~~~~~~~~~~~~~


Rayleigh-Taylor Instability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rayleigh–Taylor Instability consists of a dense fluid over a lighter fluid in a gravitational field, and is perturbed to initiate the instability.

- **Gravitational acceleration** acts in the x-direction: ``g = g_x``
- **Pressure at interface**: :math:`p_0 = 2.5`
- **Velocity perturbation (in the y-direction)**:

  .. math::

     v(x, y) = v_0 \cdot \left( \frac{1 + \cos(4\pi x)}{2} \right) \left( \frac{1 + \cos(3\pi y)}{2} \right), \quad v_0 = 0.05

- **Bottom fluid (light)**:
  - Density: :math:`\rho_L = 1.0`
  - x-velocity: :math:`u_L = 0.0`

- **Top fluid (heavy)**:
  - Density: :math:`\rho_H = 2.0`
  - x-velocity: :math:`u_H = 0.0`

The pressure is initialized to maintain hydrostatic equilibrium across the fluid interface. The pressure :math:`p(x)` is given by:

.. math::

   p(x) = 
   \begin{cases}
   p_0 + \rho_L g x, & \text{for } x < 0.75 \\
   p_0 + \rho_L g \cdot 0.75 + \rho_H g (x - 0.75), & \text{for } x \geq 0.75
   \end{cases}


At each point in the domain, the following quantities are initialized:

- **Density**:

  .. math:: 
     \rho(x) = 
     \begin{cases}
     \rho_L, & x < 0.75 \\
     \rho_H, & x \geq 0.75
     \end{cases}

- **x-Momentum density**:

  .. math:: 
     (\rho u)(x) = \rho(x) \cdot u(x) = 0

- **Gravitational source term (x-direction)**:

  .. math:: 
     (\rho g)(x) = \rho(x) \cdot g

- **Total energy density**:

  .. math::

     E(x) = \rho \cdot e(\rho, p) + \frac{1}{2} \rho \left( u^2 + v^2 \right)

  where :math:`e(\rho, p)` is the specific internal energy computed from the equation of state.

- **Radiation energy density**:

  .. math:: 
     E_{\text{rad}} = 0

The setup supports 1D, 2D, and 3D simulations. A cosine-modulated perturbation in the y-velocity initiates the RTI at the fluid interface. The pressure profile ensures that the fluid is initially in hydrostatic equilibrium.


Kelvin-Helmholtz Instability
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. vim: set tabstop=2 shiftwidth=2 expandtab fo=cqt tw=72 :
