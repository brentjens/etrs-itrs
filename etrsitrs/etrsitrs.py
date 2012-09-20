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

from numpy import array, dot, pi


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
    ParameterSet(translate_m = array([ 0.01,  0.02,  0.03]), term_d = 3.1400e-09, rotate_rad = array([-0.1, -0.2, -0.3]))
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


    ``ParameterSet`` also supports multiplication by a number and addition
    of two ``ParameterSet`` s

    >>> mm  = 0.001
    >>> mas = pi/(180.0*3600.0*1000.0)
    >>> parameters = ParameterSet(array([52.1, 49.3, -58.5])*mm,
    ...                           1.34e-9, 
    ...                           array([0.891, 5.390, -8.712])*mas)
    >>> parameters
    ParameterSet(translate_m = array([ 0.0521,  0.0493, -0.0585]), term_d = 1.3400e-09, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08]))
    >>> rates      = ParameterSet(array([0.1, 0.1, -1.8])*mm,
    ...                           0.08e-9, 
    ...                           array([0.081, 0.490, -0.792])*mas)
    >>> rates
    ParameterSet(translate_m = array([ 0.0001,  0.0001, -0.0018]), term_d = 8.0000e-11, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09]))
    >>> ref_epoch = 2000.0
    >>> parameters + rates*(2010.0 - ref_epoch)
    ParameterSet(translate_m = array([ 0.0531,  0.0503, -0.0765]), term_d = 2.1400e-09, rotate_rad = array([  8.24668072e-09,   4.98873278e-08,  -8.06342114e-08]))
    


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
        return ('ParameterSet(translate_m = %r, term_d = %.4e, rotate_rad = %r)' %
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


    def __mul__(self, number):
        return ParameterSet(translate_m = self.translate_m * number,
                            term_d      = self.term_d      * number,
                            rotate_rad  = self.rotate_rad  * number)

    def __add__(self, parameter_set):
        other = parameter_set
        return ParameterSet(translate_m = self.translate_m + other.translate_m,
                            term_d      = self.term_d      + other.term_d,
                            rotate_rad  = self.rotate_rad  + other.rotate_rad)
        



class DatumTransformation(object):
    r'''
    A datum transformation is used to transform coordinates from
    reference frame A to reference frame B. It is defined by a frame
    *from* which to transform, a frame *to* which to transform, the
    transformation parameters at the reference epoch, and their rates
    of change.

    **Parameters**
    
    from_frame : string
        The reference frame *from* which the parameters transform, for
        example 'ITRF2008'

    to_frame : string
        The reference frame *to* which the parameters transform, for
        example 'ETRF2000'

    parameters : ParameterSet
        The values of the transform parameters :math:`Tn`, :math:`D`,
        and :math:`Rn`.

    rates : ParameterSet
        The annual rates of change for the parameters :math:`Tn`,
        :math:`D`, and :math:`Rn`.
        
    ref_epoch : float
        The year to which the parameters are referenced. The
        parameters at ``epoch`` are ``parameters`` + ``rates`` *
        (epoch - ref_epoch)

    **Examples**

    >>> mm  = 0.001
    >>> mas = pi/(180.0*3600.0*1000.0)
    >>> DatumTransformation(
    ...     from_frame = 'ITRF2008', to_frame = 'ETRF2000',
    ...     parameters = ParameterSet(array([52.1, 49.3, -58.5])*mm,
    ...                               1.34e-9, 
    ...                               array([0.891, 5.390, -8.712])*mas),
    ...     rates      = ParameterSet(array([0.1, 0.1, -1.8])*mm,
    ...                               0.08e-9, 
    ...                               array([0.081, 0.490, -0.792])*mas),
    ...     ref_epoch  = 2000.0)
    DatumTransformation(from_frame = 'ITRF2008', to_frame = 'ETRF2000',
            parameters = ParameterSet(translate_m = array([ 0.0521,  0.0493, -0.0585]), term_d = 1.3400e-09, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08])),
            rates      = ParameterSet(translate_m = array([ 0.0001,  0.0001, -0.0018]), term_d = 8.0000e-11, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09])),
            ref_epoch  = 2000.0)
    '''
    def __init__(self, from_frame, to_frame, parameters, rates, ref_epoch):
        self.from_frame = from_frame
        self.to_frame   = to_frame
        self.parameters = parameters
        self.rates      = rates
        self.ref_epoch  = ref_epoch


    def __repr__(self):
        return ('''%s(from_frame = %r, to_frame = %r,
        parameters = %r,
        rates      = %r,
        ref_epoch  = %r)''' %
                (type(self).__name__, self.from_frame, self.to_frame,
                 self.parameters, self.rates, self.ref_epoch))


    def forward_transform(self, epoch, xyz):
        r'''
        '''
        pass


    def reverse_transform(self, epoch, xyz):
        r'''
        '''
        pass


    def convert_fn(self, from_frame, to_frame, epoch):
        r'''
        **Returns**
        
        A function that converts an XYZ vector from ``from_frame`` to
        ``to_frame``. If ``from_frame`` is equal to ``self.to_frame``
        and v.v., the function performs the inverse transform.
        '''
        pass



class ETRF2000(DatumTransformation):
    def __init__(self, from_frame, parameters, rates, epoch):
        self.to_frame   = 'ETRF2000'
        self.from_frame = from_frame




# Coefficients from Boucher and Altamimi (2011)
# "Memo: specifications for reference frame fixing in the analysis of
# a EUREF GPS campaign"
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
