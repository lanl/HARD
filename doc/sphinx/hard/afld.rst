.. |br| raw:: html

   <br />

.. _hard_afld:

AFLD
****

**Author:** Moon Bakaya Hazarika
**Date:** August 2025

Overview of Variables
=====================

Three new user-input variables have been implemented as part of AFLD:

1. :code:`limiter_id`: Integer varying from 0 to 4. Corresponds to the limiter relation to be used.
2. :code:`closure_id`: Integer varying from 0 to 4. Corresponds to the closure relation to be used.
3. :code:`adaptive_check`: bool-type (:code:`True` or :code:`False`).

   - If :code:`True`, the AFLD routine is used and values defined in :code:`limiter_id` and :code:`closure_id` are read.
   - If :code:`False`, the FLD routine (pre-2025 CDSS) is used.

Besides this, the **Eddington Factor** field was implemented to store the values returned by invoking the relation specified in :code:`closure_id`.

Limiter Relations
=====================

Here :math:`R` is the energy gradient ratio given as input to the function, along with the :code:`limiter_id`.
:math:`\lambda` is the limited flux value returned as output by the function.

.. math::

   R = \frac{|\nabla E_r|}{\kappa \rho E_r}

``limiter_id = 0``
------------------

No limiter applied:

.. math::

   \lambda = \frac{1}{3}

``limiter_id = 1``
------------------

Approximate Levermore–Pomraning limiter:

.. math::

   \lambda = \frac{2 + R}{6 + 3R + R^2}

``limiter_id = 2``
------------------

Bruenn limiter:

.. math::

   \lambda = \frac{1}{3 + R}

``limiter_id = 3``
------------------

Larsen’s Square Root limiter:

.. math::

   \lambda = \frac{1}{\sqrt{9 + R^2}}

``limiter_id = 4``
------------------

Minerbo limiter:

.. math::

   \lambda =
   \begin{cases}
   \dfrac{2}{3 + \sqrt{9 + 12R^2}}, & R < 1.5 \\
   \dfrac{1}{1 + R + \sqrt{1 + 2R}}, & \text{otherwise}
   \end{cases}

Closure Relations
=====================

Here :math:`\lambda` is the flux value given as input to the function, along with the :code:`limiter_id` and :code:`closure_id`.
:math:`E` is the Eddington Factor returned as output.

``closure_id = 0``
------------------

.. math::

   E = \lambda

``closure_id = 1``
------------------

.. math::

   E = \frac{1}{3}

``closure_id = 2``
------------------

.. math::

   E = 1 - 2\lambda

``closure_id = 3``
------------------

Analogous to Levermore–Pomraning type, i.e. :math:`E = \lambda + (\lambda R)^2`

- ``limiter_id = 0``

.. math::

   E = \frac{1}{3}

- ``limiter_id = 1``

.. math::

   t = \frac{1}{2}\max(0, 1-3\lambda) + \sqrt{\max(0, (1-3\lambda)(1+5\lambda))}

.. math::

   E = \lambda + t^2

- ``limiter_id = 2``

.. math::

   E = 1 - 5\lambda + 9\lambda^2

- ``limiter_id = 3``

.. math::

   E = 1 + \lambda - 9\lambda^2

- ``limiter_id = 4``

.. math::

   E =
   \begin{cases}
   1 + 3\lambda - 2\sqrt{2\lambda}, & E < \frac{2}{9} \\
   \dfrac{1}{3}, & \text{otherwise}
   \end{cases}

``closure_id = 4``
------------------

Modification of Levermore–Pomraning type: :math:`E = \frac{1}{3} + \frac{2}{3}(\lambda R)^2`

- ``limiter_id = 0``

.. math::

   E = \frac{1}{3}

- ``limiter_id = 1``

.. math::

   t = \max(0, 1-3\lambda) + \sqrt{\max(0, (1-3\lambda)(1+5\lambda))}

.. math::

   E = \frac{1}{3} + \frac{t^2}{6}

- ``limiter_id = 2``

.. math::

   E = \frac{1}{3} + \frac{2}{3}(1 - 6\lambda + 9\lambda^2)

- ``limiter_id = 3``

.. math::

   E = \frac{1}{3} + \frac{2}{3}(1 - 9\lambda^2)

- ``limiter_id = 4``

.. math::

   E =
   \begin{cases}
   \dfrac{1}{3} + \dfrac{2}{3}\left(1 + 2\lambda - 2\sqrt{2\lambda}\right), & E < \frac{2}{9} \\
   \dfrac{5}{9} - \dfrac{2}{3}\lambda, & \text{otherwise}
   \end{cases}
