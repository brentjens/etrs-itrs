from numpy import dot
from etrsitrs.parameterset import ParameterSet


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

    >>> from numpy import array
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

    >>> from numpy import array
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

    >>> from numpy import array, pi
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
    >>> itrf_to_etrf = transform.convert_fn('ITRF2008', 'ETRF2000', epoch = 2005.0)
    >>> onsala_etrf2000 = itrf_to_etrf(onsala_itrf2008)
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_etrf2000))
    3370658.848, 711876.948, 5349786.770

    Not bad at all. We also have the reverse transform, from ETRF2000 to ITRF2008:

    >>> onsala_etrf2000 = array([3370658.848, 711876.948, 5349786.770])
    >>> etrf_to_itrf = transform.convert_fn('ETRF2000', 'ITRF2008', epoch = 2005.0)
    >>> onsala_itrf2008 = etrf_to_itrf(onsala_etrf2000)
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_itrf2008))
    3370658.542, 711877.138, 5349786.952
    
    For single use, one can call the *convert* method, which under the
    hood first creates a conversion function. Note that this is fairly
    wasteful in terms of cpu cycles:

    >>> print('%.3f, %.3f, %.3f' %
    ...       tuple(transform.convert(onsala_itrf2008, 'ITRF2008', 'ETRF2000', 2005.0)))
    3370658.848, 711876.948, 5349786.770
    >>> print('%.3f, %.3f, %.3f' %
    ...       tuple(transform.convert(onsala_etrf2000, from_frame = 'ETRF2000', to_frame = 'ITRF2008', epoch = 2005.0)))
    3370658.542, 711877.138, 5349786.952

    But be careful:
    
    >>> transform.convert(onsala_etrf2000, from_frame = 'ETRF2000', to_frame = 'ITRF2005', epoch = 2005.0)
    Traceback (most recent call last):
    ...
    ValueError: No transform 'ETRF2000' -> 'ITRF2005' only 'ETRF2000' <-> 'ITRF2008'.
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
            raise ValueError('No transform %r -> %r only %r <-> %r.' %
                             (from_frame, to_frame,
                              self.to_frame, self.from_frame))
        def convert_function(xyz_m):
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
        return convert_function
        

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
