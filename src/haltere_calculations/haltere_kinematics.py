import sympy


def haltere_velocity(param, verbose=False):
    """
    Calculates symbolic expressions for the velocities of right and left
    haltere end knobs as functions of time and the basic haltere parameters.

    Parameters
    ----------
    param : dict
        dictionary containing the haltere parameters where 
        param = {
            't'     : sympy.Symbol for time,
            'L'     : sympy.Symbol for haltere length (from base to end knob CM),
            'b'     : sympy.Symbol for haltere separation (widthwise),
            'beta'  : sympy.Symbol for haltere tilt angle, 
            'theta' : sympy.Function  for haltere stroke position angle,
        }

    **kwargs (dict)
        verbose : bool
            Specifies whether or not to print verbose information (default=False)

    Returns
    -------
    vel_knob_r : sympy.Matrix  
        Right haltere end knob velocity vector, shape=(3,1)
            

    vel_know_l : sympy.Matrix 
        Left haltere end knob velocity vector, shape=(3,1)
            
    """

    # Get haltere position vectors
    pos_knob_r, pos_knob_l = haltere_position(param, verbose=False)

    # Compute velocities of right and left haltere knobs
    t = param['t']
    vel_knob_r = sympy.diff(pos_knob_r, t)
    vel_knob_l = sympy.diff(pos_knob_l, t)

    if verbose:
        print()
        print('right haltere knob velocity vector = ')
        sympy.pprint(vel_knob_r)

        print()
        print('left haltere knob velocity vector = ')
        sympy.pprint(vel_knob_l)

    return vel_knob_r, vel_knob_l



def haltere_position(param, verbose=False):
    """
    Calculates symbolic expressions for the (vectory) positions of the right
    and left haltere end knobs as a functions of time and the basic haltere
    parameters.     

    Parameters
    ----------
    param : dict
        dictionary containing the haltere parameters where 
        param = {
            't'     : sympy.Symbol for time,
            'L'     : sympy.Symbol for haltere length (from base to end knob CM),
            'b'     : sympy.Symbol for haltere separation (widthwise),
            'beta'  : sympy.Symbol for haltere tilt angle, 
            'theta' : sympy.Function  for haltere stroke position angle,
        }

    **kwargs (dict)
        verbose : bool 
            Specifies whether or not to print verbose information (default=False)

    Returns
    -------
        pos_knob_r : sympy.Matrix 
            Right haltere end knob position vector, shape=(3,1)

        pos_knob_l : sympy.Matrix
            Left haltere end knob position vector, shape=(3,1)
    """

    t = param['t']
    L = param['L']
    b = param['b']
    theta = param['theta']
    beta = param['beta']

    # Right and left haltere knob position vectors without rotation and offset
    pos0_knob_r = sympy.Matrix([L, 0, 0])
    pos0_knob_l = sympy.Matrix([-L, 0, 0])

    # Right and left haltere offset vectors
    offset_r = sympy.Matrix([ b/2, 0, 0])
    offset_l = sympy.Matrix([-b/2, 0, 0])

    # Right and left stroke rotation matrices
    rot_stroke_r = sympy.rot_axis2(theta(t))
    rot_stroke_l = sympy.rot_ccw_axis2(theta(t))

    # Right and left tilt rotation matrices
    rot_tilt_r = sympy.rot_axis3(beta)
    rot_tilt_l = sympy.rot_ccw_axis3(beta)

    # Right and left haltere knob position vectors with rotation and offset
    pos_knob_r = rot_tilt_r*rot_stroke_r*pos0_knob_r + offset_r
    pos_knob_l = rot_tilt_l*rot_stroke_l*pos0_knob_l + offset_l

    if verbose:
        print()
        print(f'stroke position angle  = {sympy.pretty(theta)}')
        print(f'caudad tilt angle      = {sympy.pretty(beta)}')

        print()
        print('right knob position vector (w/o rotation or offset) = ')
        sympy.pprint(pos0_knob_r)

        print()
        print('left knob position vector (w/o rotation or offset) = ')
        sympy.pprint(pos0_knob_l)

        print()
        print('right offset vector (unrotated) = ')
        sympy.pprint(offset_r)

        print()
        print('left offset vector (unrotated) = ')
        sympy.pprint(offset_l)

        print()
        print('right stroke rotation matrix = ')
        sympy.pprint(rot_stroke_r)

        print()
        print('left stroke rotation matrix = ')
        sympy.pprint(rot_stroke_l)

        print()
        print('right tilt rotation matrix = ')
        sympy.pprint(rot_tilt_r)

        print()
        print('left tilt rotation matrix = ')
        sympy.pprint(rot_tilt_l)

        print()
        print('right knob position vector with rotation and offset =')
        sympy.pprint(pos_knob_r)

        print()
        print('left knob position vector with rotation and offset =')
        sympy.pprint(pos_knob_l)

    return pos_knob_r, pos_knob_l


