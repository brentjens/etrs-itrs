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


'''

from numpy import array, dot


class ParameterSet(object):
    r'''
    A ParameterSet holds either the parameters :math:`Tn`, :math:`D`,
    and :math:`Rn`, or their derivatives with respect to time, in
    units of meters for :math:`Tn`, and radians for :math:`Rn`.

    **Parameters**

    translate_m : sequence of floats
        A vector containing the translation parameters T1, T2, and T3
        in units of meters.

    term_d : float
        Term D from Boucher and Altamimi. It is the cosine of
        a tiny rotation, minus 1.

    rotate_rad : sequence of floats
        The rotation parameters R1, R2, and R3 in units of radians.

    **Examples**

    >>> ParameterSet((0.01, 0.02, 0.03), 3.14e-9, [-0.1, -0.2, -0.3])
    ParameterSet(translate_m = array([ 0.01,  0.02,  0.03]), term_d = 3.14e-09, rotate_rad = array([-0.1, -0.2, -0.3]))
    >>> ParameterSet((0.01, 0.02), 3.14e-9, [-0.1, -0.2, -0.3])
    Traceback (most recent call last):
    ...
    ValueError: translate_m((0.01, 0.02)) must be e sequence of 3 floats.
    >>> ParameterSet((0.01, 0.02, 0.03), 3.14, [-0.1, -0.2, -0.3])
    Traceback (most recent call last):
    ...
    ValueError: term_d(3.14) must be a very small number.
    >>> ParameterSet((0.01, 0.02, 0.03), 3.14e-9, 15.0)
    Traceback (most recent call last):
    ...
    TypeError: object of type 'float' has no len()
    >>> ParameterSet((0.01, 0.02, 0.03), 3.14e-9, (15.0,14,13,12))
    Traceback (most recent call last):
    ...
    ValueError: rotate_rad((15.0, 14, 13, 12)) must be e sequence of 3 floats.

    '''
    def __init__(self, translate_m, term_d, rotate_rad):
        self.translate_m = array(translate_m)
        self.term_d      = term_d
        self.rotate_rad  = array(rotate_rad)

        if len(translate_m) != 3:
            raise ValueError('translate_m(%r) must be e sequence of 3 floats.' %
                             (translate_m,))
        if term_d > 1e-6:
            raise ValueError('term_d(%r) must be a very small number.' %
                             (term_d,))
        if len(rotate_rad) != 3:
            raise ValueError('rotate_rad(%r) must be e sequence of 3 floats.' %
                             (rotate_rad,))


    def __repr__(self):
        return ('ParameterSet(translate_m = %r, term_d = %r, rotate_rad = %r)' %
                (self.translate_m, self.term_d, self.rotate_rad))


    def matrix(self):
        r'''
        **Returns**

        The matrix

        .. math::

           \left(\begin{array}{ccc}
           D  & -R3 &  R2 \\
           R3 &  D  & -R1 \\ 
           -R2 &  R1 &  D
           \end{array}\right)
         
        **Examples**

        >>> ps = ParameterSet((0.01, 0.02, 0.03), 3.14e-9, [-0.1, -0.2, -0.3])
        >>> ps.matrix()
        array([[  3.14000000e-09,   3.00000000e-01,  -2.00000000e-01],
               [ -3.00000000e-01,   3.14000000e-09,   1.00000000e-01],
               [  2.00000000e-01,  -1.00000000e-01,   3.14000000e-09]])
        '''
        r_1, r_2, r_3 = self.rotate_rad
        return array([[self.term_d, -r_3       ,  r_2],
                      [r_3        , self.term_d, -r_1],
                      [-r_2       , r_1        , self.term_d]])




class CoordinateTransform(object):
    def __init__(self, from_frame, to_frame, parameters, rates, epoch):
        pass

    def __repr__():
        return 'CoordinateTransform'


class ETRF2000(CoordinateTransform):
    def __init__(self, from_frame, parameters, rates, epoch):
        self.to_frame   = 'ETRF2000'
        self.from_frame = from_frame




# Coefficients from Boucher and Altamimi (2011)
# "Memo: specifications for reference frame fixing in the analysis of a EUREF GPS campaign"
#                          |'T1' |'T2' |'T3'  |'D'  |'R1'  |'R2'  |'R3'   |
#                          |(mm) |(mm) |(mm)  |x1e-9|(mas) |(mas) |(mas)  |
PARAMETER_TABLE = [        
                   ETRF2000('ITRF2008' ,
                            [52.1, 49.3, -58.5, 1.34, 0.891, 5.390, -8.712],
                            [ 0.1,  0.1,  -1.8, 0.08, 0.081, 0.490, -0.792],
                            2000.0),
                   ETRF2000('ITRF2005' ,
                            [54.1, 50.2, -53.8, 0.40, 0.891, 5.390, -8.712],
                            [-0.2,  0.1,  -1.8, 0.08, 0.081, 0.490, -0.792],
                            2000.0),
                   ETRF2000('ITRF2000' ,
                            [54.0, 51.0, -48.0, 0.00, 0.891, 5.390, -8.712],
                            [ 0.0,  0.0,   0.0, 0.00, 0.081, 0.490, -0.792],
                            2000.0)]
