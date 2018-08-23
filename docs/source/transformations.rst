Transformations between 3D and 1D systems
=========================================

When we construct the initial conditions for a cluster simulation, we
have the choice of having a 3D or a 1D representation, if the
distribution function is spherically symmetric. Note that this does not
exclude anistropy; but if the distribution function is only a function
of energy :math:`E`, then the system is also isotropic.

The 3D representation is familiar. Position and velocity are given by
six numbers in Cartesian coordinates

.. math::
   \vec{r} = (x, y, z)\,,\quad \vec{v} = (v_x, v_y, v_z)\,.


In 1D representation the position is given by a single number

.. math::
   r = |\vec{r}| = \sqrt{x^2 + y^2 + z^2}


whereas the velocity has two components, radial and transverse

.. math::
   v_r = \vec{v}\cdot\vec{r} / \|\vec{v}\|\,,\quad
   v_t = \sqrt{\vec{v}\cdot\vec{v} - v_r^2}\,.


One way to visualize this is to think about concentric spherical shells
that contract and expand. 

We can also obtain Cartesian coordinates and velocities from :math:`r`, :math:`v_r`
and :math:`v_t`. First, we have to choose a position for our star. This
involves randomizing two angles, the polar angle :math:`\theta` and the
azimuthal angle :math:`\varphi`. The volume element in 3D spherical
coordinates imply that :math:`\varphi` should be distributed uniformly in the
interval :math:`[0, 2\pi)` (or :math:`[-\pi, pi)` or something like that, it does
not make a difference) and :math:`\cos\theta` should be distributed in
:math:`[-1, 1]` interval. We generate two independently distributed random
numbers uniform in the unit interval :math:`X_1` and :math:`X_2` and obtain

.. math::
   \phi = 2\pi X_1\,,
   \qquad \cos\theta = 2X_2-1\,,
   \quad \sin\theta = \sqrt{1-\cos^2\theta}

For the velocities, start by assuming the star is on the :math:`z`-axis.
Radial component of the velocity will be along :math:`z`-axis, and the
transverse component will be in a random direction perpendicular to
:math:`z`-axis.

.. math::
   v_z' = \pm v_r\,, v_x' = v_t \sin\alpha\,, v_y' = v_t \cos\alpha\,.

where :math:`\alpha` is a random angle in the interval :math:`[0, 2\pi)` and the
sign of :math:`v_z` is determined randomly. We then rotate this vector by angle
:math:`\theta` around :math:`y`-axis and then by angle :math:`\varphi` around :math:`z`-axis. This gives

.. math::
   \begin{split}
   \left(\begin{matrix} v_x\\v_y\\v_z\end{matrix}\right) &= 
   \left(\begin{matrix} \cos\varphi & -\sin\varphi & 0 \\
                        \sin\varphi & \cos\varphi & 0 \\
                        0 & 0 & 1\end{matrix}\right)
   \left(\begin{matrix} \cos\theta & 0 &  \sin\theta \\
                        0 & 1 & 0 \\
                        -\sin\theta & 0 & \cos\theta  \end{matrix}\right)
   \left(\begin{matrix} v_x' \\ v_y' \\ v_z'\end{matrix}\right)\\
   & = \left(\begin{matrix} 
               \sin\varphi[\cos\theta\cos\alpha \,v_t \pm \sin\theta \,v_r]
                       - \sin\varphi \sin\alpha \,v_t \\
               \cos\varphi[\cos\theta\cos\alpha \,v_t \pm \sin\theta \,v_r]
                       + \cos\varphi \sin\alpha \,v_t \\
               -\sin\theta\cos\alpha\,v_t \pm \cos\theta\,v_r
             \end{matrix}\right)
   \end{split}

These formulas can be verified by calculating
:math:`\vec{v}\cdot\vec{v} = v_t^2 + v_r^2` and
:math:`\vec{v}\cdot\vec{r} = v_r r`.
