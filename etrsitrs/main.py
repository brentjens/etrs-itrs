r'''
This module is home to the most important functions: *convert()* and
*convert_fn()*. Use *convert()* if you are only interested in one or
two coordinate conversions, but create a conversion function with the
help of *convert_fn()* if you want to convert more coordinates.

This module contains a hardcoded table of predefined
DatumTransformations, called *TRANSFORM_TABLE*. The function
*find_transform()* searches this table.

'''

import numpy
try:
    from etrsitrs.datumtransformation import DatumTransformation
    from etrsitrs.parameterset import ParameterSet
except ImportError:
    from datumtransformation import DatumTransformation
    from parameterset import ParameterSet
    



def convert_fn(from_frame, to_frame, epoch):
    r'''
    Returns a function *convert_function(xyz_m)* that converts an XYZ
    vector from ``from_frame`` to ``to_frame`` at the given ``epoch``.
    
    **Parameters**
    
    from_frame : string
        Frame from which to ransform, e.g. 'ITRF2008'.
    
    to_frame : string
        Frame to which to transform, e.g. 'ETRF2000'.
           
    epoch : number
        Epoch at which the coordinates were observed, or are
        required, in years. Example: 2013.5. 
        
    **Raises**

    KeyError
        If no appropriate transform is found

    **Returns**

    A function *f(xyz_m)* that returns a *numpy.array* of length 3.

    **Examples**
    
    >>> onsala_itrf2008 = numpy.array([3370658.542, 711877.138, 5349786.952])
    >>> fn = convert_fn('ITRF2008', 'ETRF2000', 2005.0)
    >>> print('%.3f, %.3f, %.3f' % tuple(fn(onsala_itrf2008)))
    3370658.848, 711876.948, 5349786.770

    '''
    transform = find_transform(from_frame, to_frame)
    return transform.convert_fn(from_frame, to_frame, epoch)


def convert(xyz_m, from_frame, to_frame, epoch):
    r'''
    Converts the *xyz_m* vector from *from_frame* to *to_frame*. Use
    only if one has to convert one or two coordinates. Create a
    conversion function with *convert_fn()* if you have to convert a
    large number of coordinates.

    **Parameters**
        
    xyz_m : sequence of 3 floats
        The coordinates to convert.

    from_frame : string
        Frame from which to ransform, e.g. 'ITRF2008'.

    to_frame : string
        Frame to which to transform, e.g. 'ETRF2000'.
           
    epoch : number
        Epoch at which the coordinates were observed, or are
        required, in years. Example: 2013.5. 
        
    **Raises**

    KeyError
        If no appropriate transform is found.

    **Returns**
        
    A *numpy.array* of length 3.

    **Examples**
    
    >>> onsala_itrf2008 = numpy.array([3370658.542, 711877.138, 5349786.952])
    >>> onsala_etrf2000 = convert(onsala_itrf2008, 'ITRF2008', 'ETRF2000', 2005.0)
    >>> print('%.3f, %.3f, %.3f' % tuple(onsala_etrf2000))
    3370658.848, 711876.948, 5349786.770
    '''
    transform = find_transform(from_frame, to_frame)
    return transform.convert(xyz_m, from_frame, to_frame, epoch)


class ETRF2000(DatumTransformation):
    r'''
    ETRF2000 is a subclass of *DatumTransformation* that makes it
    possible to specify the 14 parameters from Boucher and Altamimi
    (2011) in the units used in their memo. That is, first the three
    translations in mm, then term D in units of :math:`10^{-9}`
    followed by the three rotations in mas. The *to_frame* is set to
    'ETRF2000'.

    The same order and units are used for the rates.

    **Parameters**
    
    from_frame : string
        The frame from which the parameters transform,
        e.g. 'ITRF2008'.

    parameters : sequence of 7 floats
        The parameters [T1 (mm), T2 (mm), T3 (mm), D (1e-9), R1 (mas),
        R2 (mas), R3 (mas)].

    rates: sequence of 7 floats
        The annual rates of change for the parameters [T1 (mm), T2
        (mm), T3 (mm), D (1e-9), R1 (mas), R2 (mas), R3 (mas)].

    ref_epoch : float
        The reference epoch at which *parameters* are valid,
        e.g. 2000.0

    **Examples**
    
    >>> #       |'T1' |'T2' |'T3'  |'D'  |'R1'  |'R2'  |'R3'   |
    >>> #       |(mm) |(mm) |(mm)  |x1e-9|(mas) |(mas) |(mas)  | 
    >>> ETRF2000('ITRF2008' ,
    ...          [52.1, 49.3, -58.5, 1.34, 0.891, 5.390, -8.712],
    ...          [ 0.1,  0.1,  -1.8, 0.08, 0.081, 0.490, -0.792],
    ...          2000.0)
    ETRF2000(from_frame = 'ITRF2008', to_frame = 'ETRF2000',
            parameters = ParameterSet(translate_m = array([ 0.0521,  0.0493, -0.0585]), term_d = 1.3400e-09, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08])),
            rates      = ParameterSet(translate_m = array([ 0.0001,  0.0001, -0.0018]), term_d = 8.0000e-11, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09])),
            ref_epoch  = 2000.0)
   
    '''
    def __init__(self, from_frame, parameters, rates, ref_epoch):
        mm  = 0.001
        mas = numpy.pi / (180.0 * 3600.0 * 1000.0)
        super(ETRF2000, self).__init__(
            from_frame = from_frame, to_frame = 'ETRF2000',
            parameters = ParameterSet(numpy.array(parameters[0:3])*mm,
                                      parameters[3]*1e-9, 
                                      numpy.array(parameters[4:])*mas),
            rates      = ParameterSet(numpy.array(rates[0:3])*mm,
                                      rates[3]*1e-9, 
                                      numpy.array(rates[4:])*mas),
            ref_epoch  = ref_epoch)
        




# Coefficients from Boucher and Altamimi (2011)
# "Memo: specifications for reference frame fixing in the analysis of
# a EUREF GPS campaign"
#                          |'T1' |'T2' |'T3'  |'D'   |'R1'  |'R2'  |'R3'   |
#                          |(mm) |(mm) |(mm)  |x 1e-9|(mas) |(mas) |(mas)  |
TRANSFORM_TABLE = [        
                   ETRF2000('ITRF2008' ,
                            [52.1, 49.3, -58.5,  1.34, 0.891, 5.390, -8.712],
                            [ 0.1,  0.1,  -1.8,  0.08, 0.081, 0.490, -0.792],
                            2000.0),
                   ETRF2000('ITRF2005' ,
                            [54.1, 50.2, -53.8,  0.40, 0.891, 5.390, -8.712],
                            [-0.2,  0.1,  -1.8,  0.08, 0.081, 0.490, -0.792],
                            2000.0),
                   ETRF2000('ITRF2000' ,
                            [54.0, 51.0, -48.0,  0.00, 0.891, 5.390, -8.712],
                            [ 0.0,  0.0,   0.0,  0.00, 0.081, 0.490, -0.792],
                            2000.0),
                   ETRF2000('ITRF97',
                            [47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF96',
                            [47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF94',
                            [47.3, 46.7, -25.3, -1.58, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF93',
                            [76.1, 46.9, -19.9, -2.07, 2.601, 6.870, -8.412],
                            [ 2.9,  0.2,   0.6, -0.01, 0.191, 0.680, -0.862],
                            2000.0),
                   ETRF2000('ITRF92',
                            [39.3, 44.7, -17.3, -0.87, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF91',
                            [27.3, 30.7, -11.3, -2.27, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF90',
                            [29.3, 34.7,   4.7, -2.57, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0),
                   ETRF2000('ITRF89',
                            [24.3, 10.7,  42.7, -5.97, 0.891, 5.390, -8.772],
                            [ 0.0,  0.6,   1.4, -0.01, 0.081, 0.490, -0.812],
                            2000.0)]




                            


def find_transform(from_frame, to_frame):
    r'''
    Finds the appropriate *DatumTransformation* in *TRANSFORM_TABLE*.
    
    **Parameters**

    from_frame : string
        Frame from which to transform, e.g. 'ITRF2008'.

    to_frame : string
        Frame to which to transform, e.g. 'ETRF2000'

    **Returns**

    A *DatumTransformation* instance that can perform the requested
    conversion.

    **Raises**

    KeyError
        If no appropriate transform is found.

    **Examples**

    >>> find_transform('ETRF2000', 'ITRF2000')
    ETRF2000(from_frame = 'ITRF2000', to_frame = 'ETRF2000',
            parameters = ParameterSet(translate_m = array([ 0.054,  0.051, -0.048]), term_d = 0.0000e+00, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08])),
            rates      = ParameterSet(translate_m = array([ 0.,  0.,  0.]), term_d = 0.0000e+00, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09])),
            ref_epoch  = 2000.0)
    >>> find_transform('ITRF2008', 'ETRF2000')
    ETRF2000(from_frame = 'ITRF2008', to_frame = 'ETRF2000',
            parameters = ParameterSet(translate_m = array([ 0.0521,  0.0493, -0.0585]), term_d = 1.3400e-09, rotate_rad = array([  4.31968990e-09,   2.61314574e-08,  -4.22369679e-08])),
            rates      = ParameterSet(translate_m = array([ 0.0001,  0.0001, -0.0018]), term_d = 8.0000e-11, rotate_rad = array([  3.92699082e-10,   2.37558704e-09,  -3.83972435e-09])),
            ref_epoch  = 2000.0)
    >>> find_transform('ITRF1833', 'ETRF2000')
    Traceback (most recent call last):
    ...
    KeyError: "No 'ITRF1833' -> 'ETRF2000' in etrsitrs.main.TRANSFORM_TABLE."

    '''
    table  = TRANSFORM_TABLE
    frames = [to_frame, from_frame]
    
    result = [transform for transform in table
              if transform.to_frame in frames
              and transform.from_frame in frames]
    
    if len(result) == 0:
        raise KeyError('No %r -> %r in etrsitrs.main.TRANSFORM_TABLE.' %
                       (from_frame, to_frame))
    else:
        return result[0]
