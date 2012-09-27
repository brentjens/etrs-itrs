from numpy import array, dot, pi

from etrsitrs.parameterset import ParameterSet
from etrsitrs.datumtransformation import DatumTransformation



def convert_fn(from_frame, to_frame, epoch):
    r'''
    '''
    pass


def convert(from_frame, to_frame, epoch):
    r'''
    '''
    pass


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
        mas = pi / (180.0 * 3600.0 * 1000.0)
        super(ETRF2000, self).__init__(
            from_frame = from_frame, to_frame = 'ETRF2000',
            parameters = ParameterSet(array(parameters[0:3])*mm,
                                      parameters[3]*1e-9, 
                                      array(parameters[4:])*mas),
            rates      = ParameterSet(array(rates[0:3])*mm,
                                      rates[3]*1e-9, 
                                      array(rates[4:])*mas),
            ref_epoch  = ref_epoch)
        




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


def find_transform(from_frame, to_frame):
    r'''
    '''
    
