from numpy import array


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

    >>> from math import pi
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
        return ('ParameterSet(translate_m = %r, term_d = %.4e, rotate_rad = %r)'
                % (self.translate_m, self.term_d, self.rotate_rad))


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
        


