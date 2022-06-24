r'''
The ``etrsitrs`` Python module converts ITRS coordinates in various
reference frames to ETRS89 coordinates in the ETRF2000 reference frame
and vice versa. The conversions are described by 14 parameter
transforms, consisting of seven parameters and their associated rates
of change per year.

The transform and the coefficients are defined in the EUREF memo
/Specifications for reference frame fixing in the analysis of a EUREF
GPS campaign/ by Claude Boucher and Zuheir Altamimi. This module
uses version 8 bis of this memo, published on 2011-May-18.

The seven parameter transform of coordinates from frame A to frame B is
defined as

.. math::

   \left(\begin{array}{c} x_B \\ y_B \\ z_B\end{array}\right) = 
   \left(\begin{array}{c} x_A \\ y_A \\ z_A\end{array}\right) +
   \left(\begin{array}{c} T1 \\ T2 \\ T3\end{array}\right) +
   \left(\begin{array}{ccc}
       D  & -R3 &  R2 \\
       R3 &  D  & -R1 \\ 
      -R2 &  R1 &  D
   \end{array}\right)
   \left(\begin{array}{c} x_A \\ y_A \\ z_A\end{array}\right)

Table XX lists the parameters to transform from ITRFyyyy to ETRF2000.

The inverse transform, given the same parameters, is

.. math::

   \left(\begin{array}{c} x_A \\ y_A \\ z_A\end{array}\right) = 
   \left(\begin{array}{c} x_B \\ y_B \\ z_B\end{array}\right) -
   \left(\begin{array}{c} T1 \\ T2 \\ T3\end{array}\right) -
   \left(\begin{array}{ccc}
       D  & -R3 &  R2 \\
       R3 &  D  & -R1 \\ 
      -R2 &  R1 &  D
   \end{array}\right)
   \left[\left(\begin{array}{c} x_B \\ y_B \\ z_B\end{array}\right) - 
         \left(\begin{array}{c} T1 \\ T2 \\ T3\end{array}\right)\right].

Here, we used that the matrix :math:`I + M`, with :math:`I` the identity matrix,
and :math:`M` the matrix in the equations above, is a rotation
matrix. Rotation matrices are unitary, hence their inverse is equal to
their transpose.

The seven parameters :math:`Tn`, :math:`D`, and :math:`Rn` are time
dependent. To correctly perform the transformation, the parameter set
:math:`P` must first be propagated to the epoch at which the ITRF
coordinates were observed, or at which the ITRF coordinates are
desired, according to

.. math::

   P(t) = P(t_0) + \dot{P}\times(t - t_0),

where :math:`t_0` is the epoch at which the parameters are valid
(2000.0 in this case), and :math:`t` is the epoch at which the ITRF
coordinates are observed or required. Both are in units of years.


Available sub modules
---------------------

parameterset
    Organises the seven parameters necessary for the transforms.

datumtransformation
    The actual forward- and reverse transform math, and handling of
    annual change of the parameters.

main
    Contains the *convert_fn()* and *convert()* functions, as well as
    a table of predefined transforms *TRANSFORM_TABLE*, for supporting
    these functions.

'''

try:
    from importlib import metadata
except ImportError: # for Python<3.8
    import importlib_metadata as metadata

__version__ = metadata.version('etrs-itrs')

try:
    from etrsitrs.main import convert, convert_fn
except ImportError:
    from main import convert, convert_fn
