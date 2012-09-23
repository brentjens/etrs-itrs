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
        




def forward_transform(xyz_m, translate_m, rotation_matrix):
    r''' Transform *xyz_m* given a translation vector and a rotation
    matrix. Only use *translate_m* and *matrix* from the
    *ParameterSet* returned by *propagate_parameters()*. Implements

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

    
    **Parameters**
    
    xyz_m : numpy.array of length 3
        The coordinates to transform in meters.
    
    translate_m : numpy.array of length 3
        Propagated (T1, T2, T3).
    
    rotation_matrix : numpy.array of shape (3, 3)
        The rotation matrix obtained by calling the
        *ParameterSet.matrix()* method on the result of
        *propagate_parameters()*.
    
    **Returns**
    
    A numpy.array of length 3 with the transformed coordinates.

    **Examples**

    At epoch 2000.0:

    >>> parameters = ParameterSet(translate_m = array([ 0.0521, 0.0493, -0.0585]),
    ...                           term_d = 1.3400e-09,
    ...                           rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08]))
    >>> onsala_itrf2008 = array([3370658.542, 711877.138, 5349786.952])
    >>> onsala_etrf2000 = forward_transform(onsala_itrf2008,
    ...                                     parameters.translate_m,
    ...                                     parameters.matrix())
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_etrf2000))
    3370658.768, 711877.023, 5349786.816
    '''        
    return xyz_m + translate_m + dot(rotation_matrix, xyz_m)



def reverse_transform(xyz_m, translate_m, rotation_matrix):
    r'''
    The opposite of *forward_transform()*. Transform xyz given a
    translation vector and a rotation matrix. Only use *translate_m*
    and *matrix* from the *ParameterSet* returned by
    *propagate_parameters()*. Implements

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

    
    **Parameters**
    
    xyz_m : numpy.array of length 3
        The coordinates to transform in meters.
    
    translate_m : numpy.array of length 3
        Propagated (T1, T2, T3).
        
    rotation_matrix : numpy.array of shape (3, 3)
        The rotation matrix obtained by calling the
        *ParameterSet.matrix()* method on the result of
        *propagate_parameters()*.

    **Returns**

    A numpy.array of length 3 with the transformed coordinates.

    **Examples**

    **Examples**

    At epoch 2000.0:

    >>> parameters = ParameterSet(translate_m = array([ 0.0521, 0.0493, -0.0585]),
    ...                           term_d = 1.3400e-09,
    ...                           rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08]))
    >>> onsala_etrf2000 = array([3370658.768, 711877.023, 5349786.816])
    >>> onsala_itrf2008 = reverse_transform(onsala_etrf2000,
    ...                                     parameters.translate_m,
    ...                                     parameters.matrix())
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_itrf2008))
    3370658.542, 711877.138, 5349786.952

    '''
    return xyz_m - translate_m - dot(rotation_matrix, xyz_m - translate_m)




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
    >>> transform = DatumTransformation(
    ...     from_frame = 'ITRF2008', to_frame = 'ETRF2000',
    ...     parameters = ParameterSet(array([52.1, 49.3, -58.5])*mm,
    ...                               1.34e-9, 
    ...                               array([0.891, 5.390, -8.712])*mas),
    ...     rates      = ParameterSet(array([0.1, 0.1, -1.8])*mm,
    ...                               0.08e-9, 
    ...                               array([0.081, 0.490, -0.792])*mas),
    ...     ref_epoch  = 2000.0)
    >>> transform
    DatumTransformation(from_frame = 'ITRF2008', to_frame = 'ETRF2000',
            parameters = ParameterSet(translate_m = array([ 0.0521,  0.0493, -0.0585]), term_d = 1.3400e-09, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08])),
            rates      = ParameterSet(translate_m = array([ 0.0001,  0.0001, -0.0018]), term_d = 8.0000e-11, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09])),
            ref_epoch  = 2000.0)

    The parameters are only valid for the reference epoch. If one
    needs to convert coordinates at any other epoch, the parameters
    must first be propagated to that epoch with the help of the rates
    of change:

    >>> epoch = 2005.0
    >>> par_2005 = transform.propagate_parameters(epoch)
    >>> par_2005
    ParameterSet(translate_m = array([ 0.0526,  0.0498, -0.0675]), term_d = 1.7400e-09, rotate_rad = array([  6.28318531e-09,   3.80093926e-08,  -6.14355897e-08]))

    
    These propagated parameters can now be used to actually transform
    coordinates from the ITRF2008 frame to ETRF2000, at the epoch
    2005.0. Here is an example for the Onsala Space Observatory, a
    EUREF class A station. According to the EUREF web site for this
    station,
    http://www.epncb.oma.be/_productsservices/coordinates/crd4station.php?station=ONSA,
    Onsala has the following coordinates:
    
    
    =========  ======= =============================   ============================   =============================
    Frame      Epoch   X                               Y                              Z  
               (y)     (m)                             (m)                            (m)  
    =========  ======= =============================   ============================   ============================= 
    ETRF2000   2005.0  :math:`3370658.847 \pm 0.001`   :math:`711876.949 \pm 0.001`   :math:`5349786.771 \pm 0.001`
    ITRF2008   2005.0  :math:`3370658.542 \pm 0.001`   :math:`711877.138 \pm 0.001`   :math:`5349786.952 \pm 0.001`
    =========  ======= =============================   ============================   =============================

    Let's see how this works out:

    >>> onsala_itrf2008 = array([3370658.542, 711877.138, 5349786.952])
    >>> onsala_etrf2000 = forward_transform(onsala_itrf2008,
    ...                                     par_2005.translate_m,
    ...                                     par_2005.matrix())
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_etrf2000))
    3370658.848, 711876.948, 5349786.770

    Not bad at all. We also have the reverse transform, from ETRF2000 to ITRF2008:

    >>> onsala_etrf2000 = array([3370658.848, 711876.948, 5349786.770])
    >>> onsala_itrf2008 = reverse_transform(onsala_etrf2000,
    ...                                     par_2005.translate_m,
    ...                                     par_2005.matrix())
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_itrf2008))
    3370658.542, 711877.138, 5349786.952
    
    
    
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


    def propagate_parameters(self, epoch):
        r'''
        Propagate the parameter set to  *epoch*. Only use parameters
        propagated with this method to *epoch* whenever you want to do
        a coordinate conversion.

        **Parameters**

        epoch : number
            The year at which one desires the parameters,
            e.g. 2013.2.

        **Returns**
        
        A ParameterSet.
        '''
        return self.parameters + self.rates*(epoch - self.ref_epoch)


    def convert_fn(self, from_frame, to_frame, epoch):
        r'''
        Returns a function *convert(xyz_m)* that converts an XYZ
        vector from ``from_frame`` to ``to_frame``. If ``from_frame``
        is equal to ``self.to_frame`` and v.v., the function performs
        the inverse transform.

        **Parameters**

        from_frame : string
            Frame from which to ransform, e.g. 'ITRF2008'.

        to_frame : string
            Frame to which to transform, e.g. 'ETRF2000'.
           
        epoch : number
            Epoch at which the coordinates were observed, or are
            required, in years. Example: 2013.5. 
        
        **Raises**

        ValueError
            if *to_frame* or *from_frame* is not in *[self.to_frame,
            self.from_frame]*. 

        **Returns**

        A function *f(xyz_m)* that returns a *numpy.array* of length
        3.
        '''
        parameters = self.propagate_parameters(epoch)
        translate_m = parameters.translate_m
        matrix      = parameters.matrix()

        if from_frame == self.from_frame and to_frame == self.to_frame:
            transform = forward_transform
        elif from_frame == self.to_frame and to_frame == self.from_frame:
            transform = reverse_transform
        else:
            raise ValueError('Cannot transform %r to %r' % (from_frame, to_frame))
        def convert(xyz_m):
            r'''
            Convert *xyz_m* to anorther datum.
            
            **Parameters**
            
            xyz_m : numpy.array of floats of length 3
                The coordinates to transform.
            
            **Returns**

            A numpy.array of floats of length 3 containing the
            transformed coordinates.
            '''
            return transform(xyz_m, translate_m, matrix)
        return convert
        

    def convert(self, xyz_m, from_frame, to_frame, epoch):
        r''' Converts the *xyz_m* vector from *from_frame* to
        *to_frame*. If *from_frame* is equal to *self.to_frame* and
        v.v., the function performs the inverse transform. Use only if
        one has to convert one or two coordinates. Create a conversion
        function with *DatumTransformation.convert_fn()* if you have
        to convert a large number of coordinates.

        **Parameters**

        from_frame : string
            Frame from which to ransform, e.g. 'ITRF2008'.

        to_frame : string
            Frame to which to transform, e.g. 'ETRF2000'.
           
        epoch : number
            Epoch at which the coordinates were observed, or are
            required, in years. Example: 2013.5. 
        
        **Raises**

        ValueError
            if *to_frame* or *from_frame* is not in *[self.to_frame,
            self.from_frame]*. 

        **Returns**
        
        A *numpy.array* of length 3.
        '''
        return self.convert_fn(from_frame, to_frame, epoch)(xyz_m)





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
