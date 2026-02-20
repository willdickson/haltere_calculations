import sympy
from .haltere_kinematics import haltere_velocity


def coriolis_force(param, verbose=False):
    """
    Calculates symbolic expressions for the forces acting on the right and left
    halteres.

    **kwargs (dict)
        verbose : bool
            Specifies whether or not to print verbose information (default=False)

    Returns:
        force_r : sympy.Matrix
            vector of symbolic expressions for force acting on the right
            haltere end knob, shape=(3,1)

        force_l : sympy.Matrix
            vector of symbolic expressions for force acting on the left
            haltere end knob, shape=(3,1)

        w_contrib_dict : dict
            force contributions due to w_i for right and left wings where
            w_contrib_dict[side][w_i] = vector of symbolic expressions for
            force acting on the side (= left or right) wing due to angular
            velocity component w_i.
    """

    m = param['m']                   # haltere end knob mass
    v = velocity()                   # haltere end knob velocity vector
    w = w_from_param(param)          # fly's angular velocity vector
    zero = sympy.Matrix([0, 0, 0])   # the zero vector

    # Cross product matrix $\left[ \omega \right]_{\time}$
    w_cross = w_cross_matrix(param)

    # List of cross product contribution matrices (cross_1, cross_2, and cross_3)
    # for separating out contribution to cross product by vector element. 
    # e.g. w x v = w_1 * cross_1 * v + w_2 * cross_2 * v + w_3 * cross_3 *v. 
    cross_list = cross_matrices()

    # Calculate coriolis force
    force = -2*m*w.cross(v)
    force = force.expand()
    force_from_w_cross = -2*m*w_cross*v

    # Check that force from w_cross is correct
    if (force - force_from_w_cross) != zero:
        print('error: force != force_from_w_cross')

    # Calculate contributions from angular velocity components w_1, w_3 and w_3
    force_w_contrib = {w_i: -2*m*w_i*cross_i*v for w_i, cross_i in zip(w, cross_list)}
    force_from_w_contrib = sum([val for _, val in force_w_contrib.items()], zero)

    # Check that force form sum of w_i*cross_i terms is correct
    if (force - force_from_w_contrib) != zero:
        print('error: force != force_from_w_contrib')

    # Substitute in explicit expressions for haltere end knob velocity to get
    # coriolis forces on right and left halteres
    v_expr_r, v_expr_l = haltere_velocity(param, verbose=False)
    
    force_r = force
    for v_i, v_i_expr in zip(v, v_expr_r):
        force_r = force_r.subs(v_i, v_i_expr)

    force_l = force
    for v_i, v_i_expr in zip(v, v_expr_l):
        force_l = force_l.subs(v_i, v_i_expr)

    # Get right haltere force contribution contribution for the w_i components
    force_w_contrib_r = {}
    for w_i, contrib_i in force_w_contrib.items():
        contrib_i_r = contrib_i 
        for v_i, v_i_expr in zip(v, v_expr_r):
            contrib_i_r = contrib_i_r.subs(v_i, v_i_expr)
        force_w_contrib_r[w_i] = contrib_i_r

    # Get left haltere force contribution contribution for the w_i components
    force_w_contrib_l = {}
    for w_i, contrib_i in force_w_contrib.items():
        contrib_i_l = contrib_i 
        for v_i, v_i_expr in zip(v, v_expr_l):
            contrib_i_l = contrib_i_l.subs(v_i, v_i_expr)
        force_w_contrib_l[w_i] = contrib_i_l

    w_contrib_dict  = {
            'right' : force_w_contrib_r, 
            'left'  : force_w_contrib_l,
            }

    if verbose:
        print()
        print(f'haltere end knob mass = {sympy.pretty(m)}')

        print()
        print(f"fly's angular velocity vector = ")
        sympy.pprint(v)

        print()
        print('fly angular velocity vector = ')
        sympy.pprint(w)

        print()
        omega = sympy.pretty(sympy.Symbol('omega'))

        print()
        print(f'coriolis force = -2 m ({omega} x v)= ')
        sympy.pprint(force)

        for w_i, contrib_i in force_w_contrib.items():
            print()
            print(f'{sympy.pretty(w_i)} force contribution')
            sympy.pprint(contrib_i)

        print()
        print('coriolis force on right haltere = ')
        sympy.pprint(force_r)

        print()
        print('coriolis force on left haltere = ')
        sympy.pprint(force_l)

        for w_i, contrib_i in force_w_contrib_r.items():
            print()
            print(f'right haltere {sympy.pretty(w_i)} force contribution')
            sympy.pprint(contrib_i)

        for w_i, contrib_i in force_w_contrib_l.items():
            print()
            print(f'left haltere {sympy.pretty(w_i)} force contribution')
            sympy.pprint(contrib_i)

    return force_r, force_l, w_contrib_dict


def normal_coriolis_force(param, verbose=False):
    """
    Calculates symbolic expressions for the component of the forces acting on
    the right and left halteres normal to the haltere rotation axes.
    """

    # Right haltere rotation axis
    ax_r = sympy.Matrix([
        sympy.sin(param['beta']), 
        sympy.cos(param['beta']), 
        0, 
        ])

    # Left haltere rotation axis
    ax_l = sympy.Matrix([
        -sympy.sin(param['beta']), 
        sympy.cos(param['beta']), 
        0, 
        ])


    force_r, force_l, w_contrib_dict = coriolis_force(param, verbose=False)
    w = w_from_param(param)

    norm_force_r = force_r.dot(ax_r)
    norm_force_r = norm_force_r.expand()
    norm_force_r = norm_force_r.simplify()

    norm_force_l = force_l.dot(ax_l)
    norm_force_l = norm_force_l.expand()
    norm_force_l = norm_force_l.simplify()

    ax_list = [ax_r, ax_l]
    side_list = ['right', 'left']

    norm_w_contrib_dict = {side: {} for side in side_list}
    for side, ax in zip(side_list, ax_list):
        for w_i in w:
            w_i_contrib = w_contrib_dict[side][w_i]
            norm_w_i_contrib = w_i_contrib.dot(ax)
            norm_w_i_contrib = norm_w_i_contrib.expand().simplify()
            norm_w_contrib_dict[side][w_i] = norm_w_i_contrib

    if verbose:
        print()
        print('normal force on right haltere')
        sympy.pprint(norm_force_r)

        print()
        print('normal force on left haltere')
        sympy.pprint(norm_force_r)

        for w_i in w:
            for side in side_list:
                print()
                print(f'normal force contribution on {side} haltere due to {sympy.pretty(w_i)} = ')
                sympy.pprint(norm_w_contrib_dict[side][w_i])
    return norm_force_r, norm_force_l, norm_w_contrib_dict

def w_from_param(param):
    """ Return angular velocity vector w from component symbols in the param """
    w_1 = param['w_1']
    w_2 = param['w_2']
    w_3 = param['w_3']
    return sympy.Matrix([w_1, w_2, w_3])


def w_cross_matrix(param):
    """ Returns cross product matrix [w]_x """
    w_1 = param['w_1']
    w_2 = param['w_2']
    w_3 = param['w_3']
    return sympy.Matrix([[0, -w_3, w_2], [w_3, 0, -w_1], [-w_2, w_1, 0]])


def velocity():
    """ Generic velocity vector """
    v_1 = sympy.Symbol('v_1')
    v_2 = sympy.Symbol('v_2')
    v_3 = sympy.Symbol('v_3')
    return sympy.Matrix([v_1, v_2, v_3])


def cross_matrices():
    """ 
    3x3 Cross product contribution matrices,  cross_i, for separating out the
    contributions from i-th vector component in cross product. Such that

        w x v =  w_1 * cross_1 * v + w_2 * cross_2 * v + w_3 * cross_3 * v

    where w = (w_1, w_2, w_3)^T  and v = (v_1, v_2, v_3)^T. 

    Another way to look at it is 

    [w]_x = w_1 * cross_1, + w_2 * cross_2 + w_3 * cross_3
    """ 
    cross_1 = sympy.Matrix([[ 0, 0, 0], [ 0, 0,-1], [ 0, 1, 0]])
    cross_2 = sympy.Matrix([[ 0, 0, 1], [ 0, 0, 0], [-1, 0, 0]])
    cross_3 = sympy.Matrix([[ 0,-1, 0], [ 1, 0, 0], [ 0, 0, 0]])
    return [cross_1, cross_2, cross_3]



